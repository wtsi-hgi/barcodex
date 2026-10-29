"""
S1 extraction: cut the barcode / the target out of its read.

Three methods, chosen independently for the barcode and the target with
--method <barcode method> <target method>:
  none   the whole read is the sequence
  pos    by position: --*_start (0-based) and --*_len (0 = to the end of the read)
  flank  the sequence between an upstream (--*_up) and a downstream (--*_down) flanking
         sequence.  Flanks are matched on the sequence only (base qualities are not used),
         with up to --flank_mismatch Hamming mismatches.  The downstream flank is searched
         after the end of the upstream flank.  An exact flank match is preferred; otherwise
         the first position within the mismatch limit is used.

Optional filters, applied after the extraction:
  max_len            hard limit set by the code (targets: 255 bp, the width of the quality mask)
  --*_min_len        reject sequences shorter than this (e.g. --target_min_len 170)
  --*_allowed_lens   keep only sequences of these lengths, e.g. "19,20,21" or "18-22"
  --barcode_temp     barcode template, e.g. NNNNATNNNNATNN: the barcode must have the template
                     length and match its fixed (non-N) bases with <= --barcode_temp_mismatch
                     mismatches

CHECK ORDER  (a rejected read pair is counted once, under the FIRST reason that hits it, so the
reject counters add up to the number of rejected pairs; see reject_reason_order)
  1 target extraction   2 barcode extraction   3 target filters   4 barcode filters
"""

from dataclasses import dataclass, field

METHODS = ("none", "pos", "flank")
TEMPLATE_BASES = set("ACGTN")


# ---- extraction settings ---------------------------------------------------

@dataclass(frozen=True)
class ExtractSpec:
    """
    Extraction settings of one read role (barcode or target).  Plain data: it is shipped to the
    parse worker processes, which build their SequenceExtractor from it.
    Parameters:
            -- role: str: 'barcode' or 'target' (prefix of the reject reasons)
            -- method: str: 'none' (whole read), 'pos' (by position) or 'flank'
            -- start: int: pos method: 0-based start
            -- length: int: pos method: length (0 = to the end of the read)
            -- require_full_len: bool: pos method: reject reads too short for start + length
            -- up: str: flank method: upstream flanking sequence
            -- down: str: flank method: downstream flanking sequence
            -- flank_mismatch: int: flank method: max Hamming mismatches per flank
            -- max_len: int: reject longer sequences (0 = off; set by the code, not an option)
            -- min_len: int: reject shorter sequences (0 = off)
            -- allowed_lens: frozenset[int]: accepted lengths (empty = any length)
            -- template: str: template of A/C/G/T/N ('' = off)
            -- template_mismatch: int: max mismatches at the fixed bases of the template
    Returns:
            -- ExtractSpec: immutable settings object
    """
    role: str
    method: str = "none"
    start: int = 0
    length: int = 0
    require_full_len: bool = True
    up: str = ""
    down: str = ""
    flank_mismatch: int = 0
    max_len: int = 0
    min_len: int = 0
    allowed_lens: frozenset = field(default_factory=frozenset)
    template: str = ""
    template_mismatch: int = 0


def parse_allowed_lens(text):
    """
    Parse a length list such as '19,20,21', '18-22' or '16,18-20'.
    Parameters:
            -- text: str or None: comma-separated lengths and/or inclusive ranges
    Returns:
            -- frozenset[int]: accepted lengths (empty set for None / '')
    """
    lengths = set()
    for item in (text or "").split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            low, high = (int(x) for x in item.split("-", 1))
            lengths.update(range(low, high + 1))
        else:
            lengths.add(int(item))
    if any(x <= 0 for x in lengths):
        raise ValueError("lengths must be positive: %r" % text)
    return frozenset(lengths)


def spec_from_args(args, role, require_full_len, max_len=0):
    """
    Build the ExtractSpec of one role from the command line (--method and --<role>_* options).
    Parameters:
            -- args: argparse.Namespace: parsed command line (args.method = [barcode, target])
            -- role: str: 'barcode' or 'target'
            -- require_full_len: bool: pos method rejects reads shorter than start + length
            -- max_len: int: reject longer sequences (0 = off)
    Returns:
            -- ExtractSpec: settings of this role
    """
    def opt(name, default=None):
        """
        Value of the option --<role>_<name> (default when the role has no such option).
        Parameters:
                -- name: str: option name without the role prefix
                -- default: object: value for options that exist for the barcode only
        Returns:
                -- object: option value
        """
        return getattr(args, "%s_%s" % (role, name), default)
    method = args.method[0 if role == "barcode" else 1] if args.method else "none"
    return ExtractSpec(role=role, method=method, start=opt("start"), length=opt("len"),
                       require_full_len=require_full_len,
                       up=(opt("up") or "").upper(), down=(opt("down") or "").upper(),
                       flank_mismatch=args.flank_mismatch, max_len=max_len,
                       min_len=opt("min_len"),
                       allowed_lens=parse_allowed_lens(opt("allowed_lens")),
                       template=(opt("temp") or "").upper(),
                       template_mismatch=opt("temp_mismatch", 0))


def validate_spec(spec):
    """
    Check that the options of one role fit its method and are valid.  Options of another method
    are an error (not silently ignored), so a wrong option is noticed before hours of computing.
    Parameters:
            -- spec: ExtractSpec: settings to check
    Returns:
            -- list[str]: error messages (empty = valid)
    """
    errors = []
    role = spec.role
    # step1 : method and the options that belong to it
    if spec.method not in METHODS:
        errors.append("%s method must be one of %s" % (role, ", ".join(METHODS)))
    if spec.method != "pos" and (spec.start or spec.length):
        errors.append("--%s_start / --%s_len need the %s method pos" % (role, role, role))
    if spec.method != "flank" and (spec.up or spec.down):
        errors.append("--%s_up / --%s_down need the %s method flank" % (role, role, role))
    if spec.method == "flank" and not (spec.up and spec.down):
        errors.append("%s method flank needs --%s_up and --%s_down" % (role, role, role))
    # step2 : numbers and flanks
    for name, value in (("--%s_start" % role, spec.start), ("--%s_len" % role, spec.length),
                        ("--flank_mismatch", spec.flank_mismatch),
                        ("--%s_min_len" % role, spec.min_len),
                        ("--%s_temp_mismatch" % role, spec.template_mismatch)):
        if value < 0:
            errors.append("%s must be >= 0" % name)
    for name, flank in (("up", spec.up), ("down", spec.down)):
        if flank and set(flank) - TEMPLATE_BASES:
            errors.append("--%s_%s may only contain A, C, G, T, N" % (role, name))
        if flank and spec.flank_mismatch >= len(flank):
            errors.append("--flank_mismatch must be smaller than the length of --%s_%s"
                          % (role, name))
    # step3 : length filters must leave something to keep
    if spec.max_len and spec.min_len > spec.max_len:
        errors.append("--%s_min_len is longer than the %d bp limit" % (role, spec.max_len))
    if spec.allowed_lens and spec.min_len and max(spec.allowed_lens) < spec.min_len:
        errors.append("--%s_allowed_lens and --%s_min_len exclude every length" % (role, role))
    # step4 : template: A/C/G/T/N only, its length must pass the length filters
    if spec.template:
        n_template = len(spec.template)
        if set(spec.template) - TEMPLATE_BASES:
            errors.append("--%s_temp may only contain A, C, G, T, N" % role)
        if "N" not in spec.template:
            errors.append("--%s_temp has no N (no variable position)" % role)
        if spec.template_mismatch >= n_template - spec.template.count("N") > 0:
            errors.append("--%s_temp_mismatch must be smaller than the number of fixed bases "
                          "of --%s_temp" % (role, role))
        if spec.allowed_lens and n_template not in spec.allowed_lens:
            errors.append("--%s_allowed_lens excludes the --%s_temp length (%d)"
                          % (role, role, n_template))
        if spec.min_len > n_template:
            errors.append("--%s_min_len is longer than --%s_temp (%d)" % (role, role, n_template))
        if spec.method == "pos" and spec.length and spec.length != n_template:
            errors.append("--%s_len (%d) differs from the --%s_temp length (%d)"
                          % (role, spec.length, role, n_template))
    return errors


def max_length(spec, read_len):
    """
    Longest sequence the extraction can return for reads of read_len bases (used to decide
    whether barcode qualities fit in one 64-bit word).
    Parameters:
            -- spec: ExtractSpec: settings of the role
            -- read_len: int: read length
    Returns:
            -- int: upper bound of the extracted length
    """
    if spec.method == "pos":
        bound = spec.length or max(0, read_len - spec.start)
    elif spec.method == "flank":
        bound = max(0, read_len - len(spec.up) - len(spec.down))
    else:
        bound = read_len
    if spec.allowed_lens:
        bound = min(bound, max(spec.allowed_lens))
    if spec.template:
        bound = min(bound, len(spec.template))
    if spec.max_len:
        bound = min(bound, spec.max_len)
    return bound


# ---- reject reasons --------------------------------------------------------

def reject_reasons(spec, stage):
    """
    Reject reasons that can occur with these settings in one check stage, in check order.
    Parameters:
            -- spec: ExtractSpec: settings of the role
            -- stage: str: 'extract' (SequenceExtractor.span) or 'filter' (SequenceExtractor.check)
    Returns:
            -- list[tuple]: [(key, description)], key = '<role>_<reason>'
    """
    role = spec.role
    reasons = []
    if stage == "extract":
        if spec.method == "flank":
            reasons.append(("upstream_not_found", "upstream flank (--%s_up) not found" % role))
            reasons.append(("downstream_not_found",
                            "downstream flank (--%s_down) not found after the upstream flank"
                            % role))
        if spec.method == "pos" and spec.require_full_len and spec.length:
            reasons.append(("too_short", "read shorter than --%s_start + --%s_len" % (role, role)))
        reasons.append(("empty", "empty %s (read too short, or flanks adjacent)" % role))
    else:
        if spec.max_len:
            reasons.append(("too_long", "%s longer than %d bp" % (role, spec.max_len)))
        if spec.min_len:
            reasons.append(("below_min_len", "%s shorter than --%s_min_len %d"
                            % (role, role, spec.min_len)))
        if spec.allowed_lens:
            reasons.append(("length_filtered", "%s length not in --%s_allowed_lens"
                            % (role, role)))
        if spec.template:
            reasons.append(("template_mismatch", "%s does not match --%s_temp %s"
                            % (role, role, spec.template)))
    return [("%s_%s" % (role, key), "read pairs rejected: " + text) for key, text in reasons]


def reject_reason_order(barcode_spec, target_spec):
    """
    All reject reasons in the order ingest_chunk checks them: target extraction, barcode
    extraction, target filters, barcode filters (the order of associate_*.py).
    Parameters:
            -- barcode_spec: ExtractSpec: barcode settings
            -- target_spec: ExtractSpec: target settings
    Returns:
            -- list[tuple]: [(key, description)]
    """
    return (reject_reasons(target_spec, "extract") + reject_reasons(barcode_spec, "extract")
            + reject_reasons(target_spec, "filter") + reject_reasons(barcode_spec, "filter"))


# ---- sequence matching -----------------------------------------------------

def within_hamming(seq, pattern, start, max_mismatch):
    """
    True when seq[start:start+len(pattern)] differs from pattern at <= max_mismatch positions.
    Parameters:
            -- seq: bytes: sequence
            -- pattern: bytes: pattern (must fit: start + len(pattern) <= len(seq))
            -- start: int: position in seq
            -- max_mismatch: int: mismatches allowed
    Returns:
            -- bool: True = match
    """
    mismatches = 0
    for i, base in enumerate(pattern):
        if seq[start + i] != base:
            mismatches += 1
            if mismatches > max_mismatch:
                return False
    return True


def find_flank(seq, flank, max_mismatch, begin=0):
    """
    Position of a flanking sequence in seq, searching from `begin`.  An exact match is taken
    first (fast C search); otherwise the first window with <= max_mismatch Hamming mismatches.
    Parameters:
            -- seq: bytes: read sequence (upper case)
            -- flank: bytes: flanking sequence (upper case)
            -- max_mismatch: int: mismatches allowed
            -- begin: int: first position searched
    Returns:
            -- int: start position of the flank, or -1 when not found
    """
    # step1 : exact match
    pos = seq.find(flank, begin)
    if pos >= 0 or max_mismatch == 0:
        return pos
    # step2 : first window within the mismatch limit
    for pos in range(begin, len(seq) - len(flank) + 1):
        if within_hamming(seq, flank, pos, max_mismatch):
            return pos
    return -1


def matches_template(seq, start, end, fixed_bases, template_len, max_mismatch):
    """
    True when seq[start:end] has the template length and differs from the template's fixed
    (non-N) bases at <= max_mismatch positions (same rule as check_barcode of sequence_utils.py).
    Parameters:
            -- seq: bytes: read sequence
            -- start / end: int: slice of the extracted sequence
            -- fixed_bases: tuple: ((position, base byte value), ...) of the non-N template bases
            -- template_len: int: template length
            -- max_mismatch: int: mismatches allowed
    Returns:
            -- bool: True = match
    """
    if end - start != template_len:
        return False
    mismatches = 0
    for pos, base in fixed_bases:
        if seq[start + pos] != base:
            mismatches += 1
            if mismatches > max_mismatch:
                return False
    return True


# ---- per-read extraction ---------------------------------------------------

class SequenceExtractor:
    """
    Locates the barcode / target inside a read (span) and applies the filters (check) according
    to an ExtractSpec.  Built inside each parse worker (it holds bound methods and is not meant
    to be pickled).
    Parameters:
            -- spec: ExtractSpec: settings of the role
    Returns:
            -- SequenceExtractor: call span(seq), then check(seq, start, end), for every read
    """

    def __init__(self, spec):
        """
        Pre-compute the byte patterns and reject keys and choose the method.
        Parameters:
                -- spec: ExtractSpec: settings of the role
        Returns:
                -- None: constructor
        """
        self.spec = spec
        self.up = spec.up.encode()
        self.down = spec.down.encode()
        self.end = spec.start + spec.length if spec.length else None
        self.max_len = spec.max_len
        self.min_len = spec.min_len
        self.allowed_lens = spec.allowed_lens
        self.template_len = len(spec.template)
        self.fixed_bases = tuple((i, ord(b)) for i, b in enumerate(spec.template) if b != "N")
        prefix = spec.role + "_"
        self.r_up = prefix + "upstream_not_found"
        self.r_down = prefix + "downstream_not_found"
        self.r_short = prefix + "too_short"
        self.r_empty = prefix + "empty"
        self.r_too_long = prefix + "too_long"
        self.r_min_len = prefix + "below_min_len"
        self.r_length = prefix + "length_filtered"
        self.r_template = prefix + "template_mismatch"
        self._span = {"none": self._span_whole, "pos": self._span_position,
                      "flank": self._span_flank}[spec.method]

    def _span_whole(self, seq):
        """
        Method none: the whole read.
        Parameters:
                -- seq: bytes: read sequence
        Returns:
                -- tuple: (0, len(seq))
        """
        return 0, len(seq)

    def _span_position(self, seq):
        """
        Method pos: fixed window [start, start + length) (length 0 = to the end of the read).
        Parameters:
                -- seq: bytes: read sequence
        Returns:
                -- tuple or str: (start, end), or a reject key
        """
        start, end = self.spec.start, self.end
        if end is None or end > len(seq):
            if end is not None and self.spec.require_full_len:
                return self.r_short
            end = len(seq)
        return start, end

    def _span_flank(self, seq):
        """
        Method flank: window between the upstream and the downstream flank.
        Parameters:
                -- seq: bytes: read sequence (upper case)
        Returns:
                -- tuple or str: (start, end), or a reject key
        """
        mismatch = self.spec.flank_mismatch
        up = find_flank(seq, self.up, mismatch)
        if up < 0:
            return self.r_up
        start = up + len(self.up)
        end = find_flank(seq, self.down, mismatch, start)
        if end < 0:
            return self.r_down
        return start, end

    def span(self, seq):
        """
        Locate the sequence in one read (check stage 'extract').
        Parameters:
                -- seq: bytes: read sequence (upper case)
        Returns:
                -- tuple or str: (start, end) slice of the read, or the reject key
        """
        span = self._span(seq)
        if span.__class__ is not str and span[1] <= span[0]:
            return self.r_empty
        return span

    def check(self, seq, start, end):
        """
        Apply the filters to an extracted sequence (check stage 'filter').
        Parameters:
                -- seq: bytes: read sequence (upper case)
                -- start / end: int: slice returned by span
        Returns:
                -- str or None: reject key, or None when the sequence passes
        """
        length = end - start
        if self.max_len and length > self.max_len:
            return self.r_too_long
        if length < self.min_len:
            return self.r_min_len
        if self.allowed_lens and length not in self.allowed_lens:
            return self.r_length
        if self.template_len and not matches_template(
                seq, start, end, self.fixed_bases, self.template_len,
                self.spec.template_mismatch):
            return self.r_template
        return None

    def locate(self, seq):
        """
        span + check in one call (for a single role; ingest_chunk calls them separately to keep
        the check order of both roles).
        Parameters:
                -- seq: bytes: read sequence (upper case)
        Returns:
                -- tuple or str: (start, end), or the reject key
        """
        span = self.span(seq)
        if span.__class__ is str:
            return span
        return self.check(seq, *span) or span
