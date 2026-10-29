"""
S4 sequence comparison: quality-aware edit distance and 2-bit sequence packing helpers.
"""

from .constants import COST_SCALE

# popcount(x): number of 1 bits (int.bit_count needs Python >= 3.10).
try:
    popcount = int.bit_count                # Python >= 3.10
except AttributeError:
    def popcount(x):
        """
        Count the 1 bits of an integer.  Fallback used only on Python < 3.10 (no int.bit_count).
        Used by S4 (cluster_group) to count mismatching 2-bit lanes.
        Parameters:
                -- x: int: non-negative integer
        Returns:
                -- int: number of bits that are set in x
        """
        return bin(x).count("1")

def band_dp(seq_a, qual_a, seq_b, qual_b, max_hard, qual_threshold, band):
    """
    Banded edit distance with quality-aware substitution cost.
      Substitution       : costs 1 (low) if either base is low quality, else COST_SCALE (hard).
      Insertion/deletion : always COST_SCALE (hard).
    Early exit when the best cell of a row is already over the limit.
    Parameters:
            -- seq_a: bytes: sequence A (ACGT)
            -- qual_a: bytes: one value per base of A (below qual_threshold = low quality)
            -- seq_b: bytes: sequence B (ACGT)
            -- qual_b: bytes: one value per base of B
            -- max_hard: int: maximum number of hard edits
            -- qual_threshold: int: a base with quality < this is low quality
            -- band: int: half-width of the band = maximum indels; |len_a - len_b| must be <= band
    Returns:
            -- int or None: combined cost hard*COST_SCALE + low, or None when it exceeds max_hard hard
               edits
    """
    # step1 : set up the DP: cost limit, band width and the first row (j insertions)
    #      |----   substep limit = largest cost that still has <= max_hard hard edits (low edits are < COST_SCALE)
    len_a, len_b = len(seq_a), len(seq_b)
    INF = 1 << 60
    limit = max_hard * COST_SCALE + COST_SCALE - 1
    width = 2 * band + 1
    prev = [INF] * (width + 1)
    for j in range(0, min(len_b, band) + 1):                # first row: j insertions
        prev[j + band] = j * COST_SCALE
    # step2 : fill the band row by row (i = position in seq_a, j = position in seq_b)
    #      |----   substep k = j - i + band is the column inside the band; cells outside the band stay INF
    for i in range(1, len_a + 1):
        cur = [INF] * (width + 1)
        base_a = seq_a[i - 1]
        qual_ai = qual_a[i - 1]
        row_min = INF
        j_lo = i - band if i > band else 0
        j_hi = i + band if i + band < len_b else len_b
        for j in range(j_lo, j_hi + 1):
            k = j - i + band                                # column inside the band
            if j == 0:
                cost = i * COST_SCALE                       # i deletions
            else:
                #      |----   substep diagonal move: match is free; substitution is low cost if either base is low quality, else hard
                if base_a == seq_b[j - 1]:
                    cost = prev[k]                          # match
                else:
                    q_min = qual_b[j - 1]
                    if qual_ai < q_min:
                        q_min = qual_ai
                    cost = prev[k] + (1 if q_min < qual_threshold else COST_SCALE)   # substitution
                #      |----   substep vertical move (gap): always a hard edit
                insert = prev[k + 1] + COST_SCALE
                if insert < cost:
                    cost = insert
                #      |----   substep horizontal move (gap in the other sequence): always a hard edit
                if k > 0:
                    delete = cur[k - 1] + COST_SCALE
                    if delete < cost:
                        cost = delete
            cur[k] = cost
            if cost < row_min:
                row_min = cost
        #      |----   substep early exit: a row whose best cell is over the limit cannot recover
        if row_min > limit:
            return None
        prev = cur
    # step3 : read the end cell (len_a, len_b) and apply the hard-edit limit
    cost = prev[len_b - len_a + band]
    return cost if cost <= limit else None


def mate_cost(seq_a, qual_a, seq_b, qual_b, max_hard, max_low, band, qual_threshold):
    """
    (hard, low) edit counts between two sequences, or None if over the limits.
      Equal length   -> fast Hamming comparison first; the banded DP is only needed when Hamming
                        fails and indels are allowed.
      Unequal length -> banded DP (only if the length gap fits in the band).
    Parameters:
            -- seq_a: bytes: sequence A
            -- qual_a: bytes: per-base quality of A
            -- seq_b: bytes: sequence B
            -- qual_b: bytes: per-base quality of B
            -- max_hard: int: maximum hard edits
            -- max_low: int: maximum low-quality edits
            -- band: int: maximum indels (0 = Hamming only)
            -- qual_threshold: int: a base with quality < this is low quality
    Returns:
            -- tuple or None: (hard, low), or None if a limit is exceeded
    """
    # step1 : reject the pair when the length difference does not fit in the band
    len_a, len_b = len(seq_a), len(seq_b)
    if abs(len_a - len_b) > band:
        return None
    # step2 : equal length: fast quality-aware Hamming comparison
    if len_a == len_b:
        if seq_a == seq_b:
            return 0, 0
        #      |----   substep count hard / low mismatches, stop as soon as a limit is exceeded
        hard = low = 0
        for x, y, u, v in zip(seq_a, seq_b, qual_a, qual_b):
            if x != y:
                if (u if u < v else v) < qual_threshold:
                    low += 1
                else:
                    hard += 1
                if hard > max_hard or low > max_low:
                    break
        else:
            return hard, low                                # Hamming result is within limits
        #      |----   substep Hamming failed: an alignment with indels is only tried when indels / hard edits are allowed
        if band == 0 or max_hard == 0:
            return None
    # step3 : banded alignment (unequal length, or Hamming failed) and the final limit check
    cost = band_dp(seq_a, qual_a, seq_b, qual_b, max_hard, qual_threshold, band)
    if cost is None:
        return None
    hard, low = divmod(cost, COST_SCALE)
    if hard > max_hard or low > max_low:
        return None
    return hard, low


# ---- sequence packing helpers -----------------------------------------

# 4 bases per byte (2 bits each) lookup, used to unpack stored sequences.
UNPACK_TABLE = [bytes("ACGT"[(b >> 6) & 3] + "ACGT"[(b >> 4) & 3] + "ACGT"[(b >> 2) & 3]
                      + "ACGT"[b & 3], "ascii") for b in range(256)]

# lane_mask() results, cached per padded length.
LANE_MASK_CACHE = {}

# ACGT -> 0123 (digits of the 2-bit packing).
BASE_TO_DIGIT = bytes.maketrans(b"ACGT", b"0123")

def lane_mask(padded_len):
    """
    Integer 0101...01 with one 1 in the low bit of every 2-bit lane (cached per length).
    AND-ing (x | x >> 1) with it leaves one bit per mismatching base.
    Parameters:
            -- padded_len: int: sequence length rounded up to a multiple of 4
    Returns:
            -- int: lane mask of 2*padded_len bits
    """
    mask = LANE_MASK_CACHE.get(padded_len)
    if mask is None:
        mask = LANE_MASK_CACHE[padded_len] = int("01" * padded_len, 2)
    return mask


def unpack_seq(packed, length):
    """
    Unpack 2-bit packed bytes (4 bases per byte, see run_groups) into an ACGT byte string.
    Parameters:
            -- packed: bytes: packed sequence (without the trailing length byte)
            -- length: int: number of real bases (drops the padding)
    Returns:
            -- bytes: ACGT sequence of the given length
    """
    return b"".join(map(UNPACK_TABLE.__getitem__, packed))[:length]
