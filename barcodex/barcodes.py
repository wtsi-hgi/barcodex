"""
Step 2 -- processing barcodes: 
Clustering barcodes by sequence similarity and abundance, then mapping every barcode to its cluster centre.

1. The parental barcodes are the ones with count >= --barcode_ratio (the "heavy" set). 
2. Every barcode is matched to its parents by pigeonhole-seed JOIN + quality-aware Hamming (or Levenshtein for indels).

SQL builders: packed barcodes, seed join, ranking, matching queries
"""

import json
import os

from .common import log, log_sub, meta_path, open_duckdb, sql_str
from .constants import COST_SCALE, MAX_QUAL_BARCODE_LEN, MIN_SEED_LEN

# ---- SQL methods -------------------------------------------------

def packed_barcode_sql(max_len, fixed_len):
    """
    To speed up the quality-aware matching, every barcode is packed into two 64-bit words (code_expr, conf_expr)
    code_expr is the barcode and conf_expr is the low-quality mask. The SQL expressions below compute these two words
        code_expr: 2 bits per base (A=0 C=1 G=2 T=3)
        conf_expr: one bit, 1 = confident, 0 = low quality

    Note: first base in the highest lane
    
    hard = bit_count(((x | x >> 1) & conf1 & conf2)),  x = code1 xor code2.

    All barcodes must have the same length (fixed_len) for the shifts to be constants (much faster).

    Parameters:
        -- max_len: int: longest barcode (at most MAX_QUAL_BARCODE_LEN = 32)
        -- fixed_len: int: common length of all barcodes, or 0 if the lengths differ
    Returns:
        -- tuple: (code_expr, conf_expr), two SQL expression strings returning UBIGINT
    """
    # step1 : one term per base position: a 2-bit code lane and a 1-bit confidence lane
    code_terms, conf_terms = [], []
    for i in range(max_len):
        #    |---- substep constant shifts when all barcodes have the same length, else shifts computed from length(b)
        if fixed_len:
            bit_shift = "%d" % (fixed_len - 1 - i)
            lane_shift = "%d" % (2 * (fixed_len - 1 - i))
        else:
            bit_shift = "greatest(length(b) - 1 - %d, 0)" % i
            lane_shift = "(2 * %s)" % bit_shift
        #    |---- substep code lane: A=0 C=1 G=2 T=3 (anything else = 0)
        code = ("(CASE substr(b, %d, 1) WHEN 'C' THEN 1 WHEN 'G' THEN 2 WHEN 'T' THEN 3 ELSE 0 END)"
                "::UBIGINT << (%s)::UBIGINT" % (i + 1, lane_shift))
        #    |---- substep confidence lane: 1 when the mask bit of the base is 0 (confident)
        conf = ("(1::UBIGINT - ((bm >> (%s)::UBIGINT) & 1::UBIGINT)) << (%s)::UBIGINT"
                % (bit_shift, lane_shift))
        #    |---- substep variable length: positions beyond the barcode end contribute 0
        if not fixed_len:
            code = "CASE WHEN length(b) > %d THEN %s ELSE 0::UBIGINT END" % (i, code)
            conf = "CASE WHEN length(b) > %d THEN %s ELSE 0::UBIGINT END" % (i, conf)
        code_terms.append("(%s)" % code)
        conf_terms.append("(%s)" % conf)
    
    # step2 : OR all lane terms together into two 64-bit words
    return (
        "(%s)::UBIGINT" % " | ".join(code_terms),
        "(%s)::UBIGINT" % " | ".join(conf_terms)
    )


def seed_join_parts(src, seed_len, n_seeds, band):
    """
    Shared logic used by both barcode-matching queries (pigeonhole seed join).

    Each barcode is divided into n_seeds non-overlapping sections (seeds), with each seed containing seed_len bases

    The pigeonhole principle guarantees that if two barcodes differ at no more than (n_seeds - 1) positions, at least 
    one corresponding seed must be identical.

    When insertions or deletions (indels) are allowed, the position of a seed may shift relative to its expected 
    position. We search for the matching seed at every possible offset within band bases of its expected position.

    For faster searching and lower memory usage, each seed is encoded as a 64-bit integer hash rather than stored 
    and compared as a DNA string. This keeps the seed index compact enough to handle up to 10^8 heavy barcodes.

    Hash collisions are possible: two different seeds may produce the same hash. This does not cause incorrect matches 
    because a hash match only generates a candidate pair. Every candidate is subsequently checked using the actual 
    barcode distance, which rejects false matches caused by collisions.

    Parameters:
        -- src: str: SQL table/subquery of the children (columns b, n, ...), aliased c
        -- seed_len: int: bases per seed
        -- n_seeds: int: number of seed slots (max edits + 1)
        -- band: int: maximum indels (0 = no indel search)
    Returns:
        -- tuple: (from_clause, seed_hash_expr, shift_filter, length_filter) SQL fragments
    """
    # with indels (band > 0): probe the child at every shift [-band .. +band] and compare lengths
    if band > 0:
        from_clause = ("%s c CROSS JOIN range(%d) AS s(slot) CROSS JOIN range(%d, %d) AS h(sh)"
                       % (src, n_seeds, -band, band + 1))
        pos = "s.slot * %d + h.sh" % seed_len
        shift_filter = "AND %s >= 0 AND %s + %d <= length(c.b)" % (pos, pos, seed_len)
        seed_expr = "hash(substr(c.b, %s + 1, %d))" % (pos, seed_len)
        length_filter = "AND abs(length(c.b) - length(p.b)) <= %d" % band
    # without indels: only the exact slot positions, no shift / length filters
    else:
        from_clause = "%s c CROSS JOIN range(%d) AS s(slot)" % (src, n_seeds)
        seed_expr = "hash(substr(c.b, s.slot * %d + 1, %d))" % (seed_len, seed_len)
        shift_filter = length_filter = ""

    return from_clause, seed_expr, shift_filter, length_filter


def rank_key_sql(cost_expr, parent_count_expr, tiebreak_expr):
    """
    Build a single HUGEINT ranking key to select the best parent barcode.

    Candidates are ranked using three criteria, in this priority order:
    1. Edit cost: prefer the parent with fewer mismatches
    2. Parent read count: if edit costs are equal, prefer parent supported by more sequencing reads.
    3. Tie-break hash: if both values are equal, prefer parent with the smaller hash

    The three values are packed into separate bit ranges of one integer:
    [ edit cost: 32 bits | inverted read count: 32 bits | hash: 64 bits ]

    The edit cost occupies the highest bits, so it always takes priority. The read count is inverted 
    (maximum count - actual count) so that parents with more reads receive smaller ranking values.
    The hash occupies the lowest bits and is only used when the first two criteria are equal.

    This allows SQL arg_min() to select the best parent using a single
    integer comparison instead of sorting by three separate columns.    

    R = C * 2^96 + (2^32 - 1 - min(N, 2^32 - 1)) * 2^64 + H
    R -- Ranking Key
    C -- edit cost
    N -- parent read count
    H -- hash of the parent
    
    Parameters:
        -- cost_expr: str: SQL expression, edit cost (smaller is better)
        -- parent_count_expr: str: SQL expression, reads of the parent (larger is better)
        -- tiebreak_expr: str: SQL expression returning a UBIGINT hash of the parent
    Returns:
        -- str: SQL expression
    """
    return (
        "(%s)::HUGEINT * 79228162514264337593543950336 + "
        "(4294967295 - LEAST(%s, 4294967295))::HUGEINT * 18446744073709551616 + "
        "(%s)::HUGEINT" % (cost_expr, parent_count_expr, tiebreak_expr)
    )


def matched_query_plain(src, seed_len, n_seeds, ratio, max_dist, band):
    """
    Find the best matching heavy parent barcode for each barcode using sequence

    Matching workflow:
    1. Generate candidate parents using a pigeonhole seed join. Each child barcode is split into seeds, 
       and matching seed hashes are used to find potentially similar heavy barcodes.
    2. Retrieve the candidate parent barcodes from the heavy table using their integer IDs (hid).
    3. Apply parent eligibility rules:
        - Parent and child must be different barcodes.
        - Parent read count must be at least X times the child read count.
        - If parent and child have equal read counts, the parent must have a lexicographically smaller barcode 
          than the child. This ensures deterministic direction of matching.
    4. Calculate the sequence distance between each candidate pair:
        - Use Hamming distance when band = 0 (substitutions only)
        - Use Levenshtein distance when band > 0 (indels allowed)
    5. Discard candidates whose sequence distance exceeds max_dist
    6. Select the best remaining parent for each child using arg_min, prioritising:
        - Smallest sequence distance
        - Largest parent read count
        - Smallest parent barcode hash as a deterministic tie-breaker

    Candidate generation uses integer seed hashes and integer parent IDs for the joins. Sequence comparisons and 
    other matching conditions are applied as filters after candidate generation.

    Parameters:
        -- src: str: SQL table/subquery with columns b, n
        -- seed_len: int: bases per seed
         -- n_seeds: int: number of seed slots
        -- ratio: float: parent count must be >= ratio x child count
        -- max_dist: int: maximum Hamming (or Levenshtein when band > 0) distance
        -- band: int: maximum indels (0 = Hamming only)
    Returns:
        -- str: SQL query returning (barcode, n, parent_id)
    """
    from_clause, seed, shift_filter, length_filter = seed_join_parts(src, seed_len, n_seeds, band)
    
    # Levenshtein when indels are allowed, else Hamming (99 = different length = rejected)
    if band > 0:
        dist = "levenshtein(c.b, p.b)"
    else:
        dist = "CASE WHEN length(c.b) = length(p.b) THEN hamming(c.b, p.b) ELSE 99 END"
    
    # ranking key for arg_min: smaller distance, then larger parent count, then hash
    key = rank_key_sql(dist, "p.n", "hash(p.b)")

    # assemble: seed join -> heavy lookup -> coverage / tie rules -> distance filter -> best parent
    return (
        "SELECT c.b AS barcode, c.n AS n, arg_min(p.hid, %s) AS parent_id FROM %s "
        "JOIN hseed q ON q.slot = s.slot AND q.seed_hash = %s "
        "JOIN heavy p ON p.hid = q.hid "
        "WHERE p.b <> c.b AND p.n >= %r * c.n "
        "AND (p.n > c.n OR (p.n = c.n AND p.b < c.b)) %s %s AND (%s) <= %d "
        "GROUP BY c.b, c.n" % (key, from_clause, seed, float(ratio),
                               shift_filter, length_filter, dist, max_dist)
    )


def matched_query_quality(src, seed_len, n_seeds, ratio, max_hard, max_low, band):
    """
    Build a SQL query for quality-aware barcode matching.

    Each barcode is matched to a more abundant parent barcode, allowing mismatches at low-quality positions 
    while applying stricter limits to mismatches at high-quality positions.

    The matching procedure has three layers:

    1. Candidate generation:
       Use seed-based joins to identify potential child-parent pairs. Count the total number of mismatches and 
       the number occurring at high-quality positions. If indels are allowed, calculate Levenshtein distance.
    2. Match evaluation:
       Convert mismatch counts into costs and reject pairs that exceed the allowed mismatch thresholds. High-quality 
       mismatches receive a bigger penalty than low-quality mismatches. A match is accepted only if both mismatch 
       limits are satisfied.
    3. Parent selection:
       Retain accepted pairs and select the best parent for each barcode using the ranking key.

    Parameters:
        -- src: str: SQL table/subquery with columns b, n, code, conf
        -- seed_len: int: bases per seed
        -- n_seeds: int: number of seed slots (max_hard + max_low + 1)
        -- ratio: float: parent count must be >= ratio x child count
        -- max_hard: int: maximum mismatches at confident bases (--barcode_dist)
        -- max_low: int: maximum extra mismatches at low-quality bases
        -- band: int: maximum indels (0 = none)
    Returns:
        -- str: SQL query returning (barcode, n, parent_id)
    """
    # Step 1: Generate seed-based candidate child-parent pairs.
    # Only pairs sharing a compatible seed are considered, avoiding an exhaustive comparison between all barcodes.
    from_clause, seed, shift_filter, length_filter = seed_join_parts(src, seed_len, n_seeds, band)

    # Step 2: Define mismatch calculations using packed barcode representations.
    # The XOR operation identifies nucleotide differences between two barcodes. Each nucleotide is encoded using 2 bits. 
    # Combining XOR result with its right-shifted copy ensures that either differing bit marks the corresponding nucleotide 
    # as a mismatch. Intersecting with both confidence masks counts only mismatches where both barcodes have high-quality bases.
    same_len = "length(c.b) = length(p.b)"
    xor_expr = "xor(c.code, p.code)"

    hard_expr = (
        "bit_count((((%s) | ((%s) >> 1::UBIGINT)) & c.conf & p.conf)::BIGINT)"
        % (xor_expr, xor_expr)
    )

    lev_col = ", levenshtein(c.b, p.b) AS lev" if band > 0 else ""

    # Step 3: Build candidate pairs and calculate their mismatch statistics.
    # Each candidate must:
    #   - Have a different barcode sequence from the child.
    #   - Meet the minimum parent-to-child read-count ratio.
    #   - Follow the abundance-based ordering used to avoid duplicate child-parent comparisons.
    layer_a = (
        "SELECT c.b AS barcode, c.n AS n, p.hid AS parent_id, p.n AS parent_n, "
        "hash(p.b) AS parent_hash, "
        "CASE WHEN %s THEN hamming(c.b, p.b) END AS total, "
        "CASE WHEN %s THEN %s END AS hard%s "
        "FROM %s JOIN hseed q ON q.slot = s.slot AND q.seed_hash = %s "
        "JOIN heavy p ON p.hid = q.hid "
        "WHERE p.b <> c.b AND p.n >= %r * c.n "
        "AND (p.n > c.n OR (p.n = c.n AND p.b < c.b)) %s %s"
        % (
            same_len, same_len, hard_expr, lev_col,
            from_clause, seed, float(ratio),
            shift_filter, length_filter
        )
    )

    # Step 4: Evaluate candidate pairs using two independent matching criteria.
    # Hamming cost:
    #   Accept only if both the high-quality mismatch count and the low-quality mismatch count are within limits.
    #   low-quality mismatches = total mismatches - high-quality mismatches
    # Levenshtein cost:
    #   Accept pairs whose edit distance is within max_hard.
    ham_cost = (
        "CASE WHEN total IS NOT NULL AND hard <= %d AND total - hard <= %d "
        "THEN hard * %d + (total - hard) END"
        % (max_hard, max_low, COST_SCALE)
    )

    lev_cost = (
        "CASE WHEN lev <= %d THEN lev * %d END" % (max_hard, COST_SCALE)
        if band > 0 else "NULL"
    )

    layer_b = (
        "SELECT barcode, n, parent_id, parent_n, parent_hash, "
        "%s AS ham_cost, %s AS lev_cost FROM (%s)"
        % (ham_cost, lev_cost, layer_a)
    )

    # Step 5: Select the best accepted parent for each child barcode.
    # A candidate is accepted if it passes either the Hamming or Levenshtein matching criterion.
    #
    # When both costs are available, use the lower cost. COALESCE assigns a large penalty to a criterion 
    # that did not accept the candidate.
    #
    # rank_key_sql adds parent abundance and a deterministic tie-breaker to the matching cost. 
    # arg_min then selects the parent with the best overall ranking for each barcode.
    best_cost = "LEAST(COALESCE(ham_cost, 2147483647), COALESCE(lev_cost, 2147483647))"

    key = rank_key_sql(best_cost, "parent_n", "parent_hash")

    return (
        "SELECT barcode, n, arg_min(parent_id, %s) AS parent_id FROM (%s) "
        "WHERE ham_cost IS NOT NULL OR lev_cost IS NOT NULL GROUP BY barcode, n"
        % (key, layer_b)
    )


# ---- DuckDB utilities -------------------------------------------------

def spill_gb(work):
    """
    Size of the DuckDB spill directory of S2 (disk used to reduce memory usage).
    Parameters:
        -- work: str: work directory (spill directory is <work>/duck_bc)
    Returns:
        -- float: size in GB
    """
    total = 0
    for d, _, files in os.walk(os.path.join(work, "duck_bc")):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(d, f))
            except OSError:
                pass
    return total / 1e9


def pin_join_order(con, enabled):
    """
    Execute the joins in the order they are written (enabled=True) or let DuckDB choose.

    Parameters:
        -- con: duckdb connection
        -- enabled: bool: True = pin the join order, False = reset to the optimizer default
    Returns:
        -- None: nothing
    """
    try:
        con.execute("SET disabled_optimizers='join_order'" if enabled else "RESET disabled_optimizers")
    except Exception as e:
        log_sub("NOTE: could not %s join-order pinning (%s)" % ("enable" if enabled else "reset", str(e).splitlines()[0]))


def dump_plan(con, query, path):
    """
    Write DuckDB's EXPLAIN of a query to a text file (to inspect the join order / build side).
    Parameters:
            -- con: duckdb connection
            -- query: str: SQL query to explain
            -- path: str: text file to write
    Returns:
            -- None: nothing (a failure only prints a note)
    """
    try:
        rows = con.execute("EXPLAIN " + query).fetchall()
        with open(path, "w") as fo:
            for row in rows:
                fo.write(str(row[-1]) + "\n")
        log_sub("query plan written to %s" % path)
    except Exception as e:
        log_sub("NOTE: could not write the query plan (%s)" % str(e).splitlines()[0])


# ---- stage driver -----------------------------------------------------

def stage_barcodes(args, work, threads, mem_gb):
    """
    S2: barcode census -> centre barcode for every barcode  (<work>/bcmap.parquet).
    Parents are barcodes with count >= --barcode_ratio (the 'heavy' set).  Every barcode is matched
    to heavy parents by seed join + quality-aware Hamming (or Levenshtein for indels); best parent =
    smallest (hard, low) distance, then largest count; chains are resolved by pointer jumping on the
    heavy set only.  Also writes <out>.meta/raw_stats.json, <out>.meta/barcode_clusters.parquet,
    the optional <out>.barcode_map.tsv and the marker <work>/bcmap.done.
    Parameters:
            -- args: argparse.Namespace: barcode_ratio/dist/max_lowq/indel, min_*_reads, max_cand,
               barcode_map, out
            -- work: str: work directory
            -- threads: int: DuckDB threads
            -- mem_gb: float: RAM budget in GB (85% is given to DuckDB)
    Returns:
            -- None: nothing; results are files
    """
    # step1 : file names, tuple source and DuckDB connection
    census_base = work + "/census_base.parquet"
    census_q = work + "/census_q.parquet"
    bcmap = work + "/bcmap.parquet"
    tuples = "read_parquet(%s)" % sql_str(work + "/tuples/t*.parquet")
    log("S2: barcode census + clustering (DuckDB, out-of-core) ...")
    log_sub("DuckDB budget for S2: %.1f GB RAM, %d thread(s); spill directory %s/duck_bc"
            % (mem_gb * 0.85, threads, work))
    #      |----   substep DuckDB with 85% of the RAM budget and its own spill directory
    con = open_duckdb(work, threads, mem_gb * 0.85, "bc", work + "/bc.duckdb")

    # step2 : census: reads per distinct barcode + AND of its low-quality masks
    con.execute("COPY (SELECT barcode AS b, SUM(cnt)::BIGINT AS n, bit_and(bm) AS bm "
                "FROM %s GROUP BY barcode) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)"
                % (tuples, sql_str(census_base)))
    #      |----   substep reads, distinct barcodes and the barcode length range
    n_distinct, n_reads, min_len, max_len = con.execute(
        "SELECT COUNT(*), SUM(n), MIN(length(b)), MAX(length(b)) FROM read_parquet(%s)"
        % sql_str(census_base)).fetchone()
    log_sub("%d reads, %d distinct barcodes (length %d-%d)" % (n_reads, n_distinct, min_len, max_len))
    #      |----   substep distinct raw target sequences (for the stats file)
    n_raw_targets = con.execute("SELECT COUNT(*) FROM (SELECT target FROM %s GROUP BY target)"
                                % tuples).fetchone()[0]
    log_sub("%d distinct raw target sequences" % n_raw_targets)

    # step3 : settings: ratio, distances, quality mode and seed layout
    ratio = max(args.barcode_ratio, 1.0)        # parent count >= ratio x child count
    max_hard = args.barcode_dist                # edits at confident bases
    band = min(args.barcode_indel, max_hard)    # indels allowed
    #      |----   substep quality-aware matching needs barcodes <= 32 bp (packed in one 64-bit word)
    use_quality = args.barcode_max_lowq > 0 and max_len <= MAX_QUAL_BARCODE_LEN
    max_low = args.barcode_max_lowq if use_quality else 0
    if args.barcode_max_lowq > 0 and not use_quality:
        log_sub("NOTE: barcodes longer than %d bp: quality-aware clustering off for barcodes"
            % MAX_QUAL_BARCODE_LEN)
    #      |----   substep pigeonhole seeds: pairs with <= max_hard+max_low differences share one of (max_hard+max_low+1)
    #      |----   substep seeds exactly; lower max_low while the seeds would get shorter than MIN_SEED_LEN
    while use_quality and max_low > 0 and min_len // (max_hard + max_low + 1) < MIN_SEED_LEN:
        max_low -= 1
        log_sub("NOTE: barcodes too short for the requested --barcode_max_lowq; using %d" % max_low)
    #      |----   substep final seed layout
    n_seeds = max_hard + max_low + 1
    seed_len = min_len // n_seeds
    if seed_len < MIN_SEED_LEN:
        raise SystemExit("ERROR: barcodes too short for --barcode_dist %d" % max_hard)
    #      |----   substep read-count thresholds for dropping light barcodes
    min_barcode_reads = (args.min_barcode_reads if args.min_barcode_reads is not None
                         else args.min_pair_reads)
    min_child_reads = max(1, args.min_child_reads)

    # step4 : census table used by all later queries
    if use_quality:
        # add code / conf words (see packed_barcode_sql); mask bit 1 = low quality
        code_expr, conf_expr = packed_barcode_sql(
            max_len, max_len if min_len == max_len else 0)
        con.execute("COPY (SELECT b, n, %s AS code, %s AS conf FROM read_parquet(%s)) "
                    "TO %s (FORMAT PARQUET, COMPRESSION ZSTD)"
                    % (code_expr, conf_expr, sql_str(census_base), sql_str(census_q)))
        census_file = census_q
        columns = "b, n, code, conf"
        log_sub("quality-aware matching: <= %d hard + <= %d low-quality mismatches, %d seeds of %d bp"
            % (max_hard, max_low, n_seeds, seed_len))
    else:
        census_file = census_base
        columns = "b, n"

    #      |----   substep helper: SQL that matches a barcode table to its best heavy parents
    def matched(src):
        """
        SQL that matches every row of `src` to its best heavy parent (quality-aware or sequence-only,
        depending on the settings of the enclosing stage_barcodes call).
        Parameters:
                -- src: str: SQL table/subquery with columns b, n (+ code, conf when quality-aware)
        Returns:
                -- str: SQL query returning (barcode, n, parent_id)
        """
        if use_quality:
            return matched_query_quality(src, seed_len, n_seeds, ratio, max_hard, max_low, band)
        return matched_query_plain(src, seed_len, n_seeds, ratio, max_hard, band)

    #      |----   substep census source for the following queries
    census = "read_parquet(%s)" % sql_str(census_file)

    # step5 : heavy barcodes = the only possible parents; build their seed index
    # Every heavy barcode gets an integer id (hid).  All big joins below use ids and 64-bit
    # seed hashes instead of the barcode strings, so the join tables stay small and in memory.
    #      |----   substep heavy table: integer id (hid) + packed barcode columns
    con.execute("CREATE OR REPLACE TABLE heavy AS SELECT ROW_NUMBER() OVER () AS hid, %s "
                "FROM %s WHERE n >= %r" % (columns, census, float(ratio)))
    n_heavy = con.execute("SELECT COUNT(*) FROM heavy").fetchone()[0]
    log_sub("%d heavy barcodes (count >= %g) can be cluster centres/parents" % (n_heavy, ratio))
    #      |----   substep seed index of the heavy set: (slot, seed hash, heavy id); seeds shared by more than --max_cand
    #      |----   substep heavy barcodes (low complexity) are dropped to bound the join size
    con.execute("""
        CREATE OR REPLACE TABLE hseed AS
        SELECT hid, slot, seed_hash FROM (
            SELECT h.hid AS hid, s.slot AS slot,
                   hash(substr(h.b, s.slot * %d + 1, %d)) AS seed_hash
            FROM heavy h CROSS JOIN range(%d) AS s(slot))
        QUALIFY count(*) OVER (PARTITION BY slot, seed_hash) <= %d
    """ % (seed_len, seed_len, n_seeds, args.max_cand))
    n_seed_rows = con.execute("SELECT COUNT(*) FROM hseed").fetchone()[0]
    log_sub("seed index built: %d entries (%d seeds of %d bp per heavy barcode); "
            "spill so far %.1f GB" % (n_seed_rows, n_seeds, seed_len, spill_gb(work)))

    # step6 : heavy barcodes: best heavy parent, then pointer jumping to the root
    #      |----   substep pin the join order: children JOIN integer seed index JOIN heavy strings
    pin_join_order(con, True)
    #      |----   substep best heavy parent of every heavy barcode (hroot: hid -> root_id)
    con.execute("CREATE OR REPLACE TABLE hroot AS "
                "SELECT h.hid AS hid, COALESCE(m.parent_id, h.hid) AS root_id "
                "FROM heavy h LEFT JOIN (%s) m ON m.barcode = h.b" % matched("heavy"))
    n_with_parent = con.execute("SELECT COUNT(*) FROM hroot WHERE root_id <> hid").fetchone()[0]
    log_sub("%d of %d heavy barcodes have a more abundant heavy parent; spill so far %.1f GB"
            % (n_with_parent, n_heavy, spill_gb(work)))
    #      |----   substep pointer jumping: root := root of root, until nothing changes
    jump_round = 0
    while True:                                 # root := root of root, until nothing changes
        jump_round += 1
        con.execute("CREATE OR REPLACE TABLE hroot2 AS "
                    "SELECT a.hid AS hid, COALESCE(b.root_id, a.root_id) AS root_id "
                    "FROM hroot a LEFT JOIN hroot b ON a.root_id = b.hid")
        n_changed = con.execute("SELECT COUNT(*) FROM hroot a JOIN hroot2 b USING (hid) "
                                "WHERE a.root_id <> b.root_id").fetchone()[0]
        con.execute("DROP TABLE hroot")
        con.execute("ALTER TABLE hroot2 RENAME TO hroot")
        log_sub("pointer jumping round %d: %d heavy barcodes changed root" % (jump_round, n_changed))
        if not n_changed:
            break
    #      |----   substep centre barcode (string) of every heavy barcode (hcentre)
    con.execute("CREATE OR REPLACE TABLE hcentre AS "
                "SELECT r.hid AS hid, r.root_id AS root_id, c.b AS centre "
                "FROM hroot r JOIN heavy c ON c.hid = r.root_id")

    # step7 : every barcode -> centre barcode (writes bcmap.parquet)
    #  heavy barcode            -> its root
    #  light barcode, parent    -> root of its heavy parent
    #  light barcode, no parent -> its own centre, only if it has >= min_barcode_reads reads
    # The expensive light-barcode join runs ONCE (table lmatch); later steps reuse it.
    # Light barcodes below --min_child_reads are ignored completely (not matched, not kept).
    #      |----   substep children = light barcodes (n < ratio and n >= --min_child_reads)
    children = "(SELECT %s FROM %s WHERE n < %r AND n >= %d)" % (columns, census, float(ratio),
                                                                  min_child_reads)
    #      |----   substep save the query plan to check the join order
    dump_plan(con, matched(children), work + "/explain_light_match.txt")
    log_sub("matching light barcodes to heavy parents (the longest step) ...")
    #      |----   substep match light barcodes to heavy parents (expensive join; run once, reused below)
    con.execute("CREATE OR REPLACE TABLE lmatch AS "
                "SELECT barcode, n, parent_id, hash(barcode) AS hk FROM (%s)" % matched(children))
    #      |----   substep release the join-order pin
    pin_join_order(con, False)
    n_lmatch = con.execute("SELECT COUNT(*) FROM lmatch").fetchone()[0]
    log_sub("%d light barcodes matched to a heavy parent; spill so far %.1f GB"
            % (n_lmatch, spill_gb(work)))
    #      |----   substep light barcodes without a parent keep their own centre when they have enough reads
    #      |----   substep (anti-join on the 64-bit barcode hash: no strings in the join table)
    con.execute("""
        CREATE OR REPLACE TABLE lown AS
        SELECT s.b AS b, s.n AS n FROM %s s
        WHERE s.n >= %d AND s.n < %r
          AND NOT EXISTS (SELECT 1 FROM lmatch m WHERE m.hk = hash(s.b))
    """ % (census, max(min_barcode_reads, min_child_reads), float(ratio)))
    n_lown = con.execute("SELECT COUNT(*) FROM lown").fetchone()[0]
    log_sub("%d light barcodes without a parent kept as their own cluster" % n_lown)
    #      |----   substep bcmap.parquet = heavy -> root | light with parent -> root of parent | light own centre
    con.execute("""
        COPY (
            SELECT h.b AS barcode, c.centre AS centre, h.n AS n         -- heavy barcodes
            FROM heavy h JOIN hcentre c ON c.hid = h.hid
            UNION ALL
            SELECT m.barcode, c.centre, m.n                              -- light, has a parent
            FROM lmatch m JOIN hcentre c ON c.hid = m.parent_id
            UNION ALL
            SELECT b, b, n FROM lown                                     -- light, own centre
        ) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)
    """ % sql_str(bcmap))
    #      |----   substep report mapped / dropped barcodes
    n_mapped, n_centres, n_mapped_reads = con.execute(
        "SELECT COUNT(*), COUNT(*) FILTER (WHERE barcode = centre), SUM(n) FROM read_parquet(%s)"
        % sql_str(bcmap)).fetchone()
    log_sub("%d distinct barcodes -> %d barcode clusters; dropped %d barcodes (%d reads) "
            "(unmatched with < %d reads, or < %d reads = --min_child_reads)"
            % (n_mapped, n_centres, n_distinct - n_mapped, n_reads - (n_mapped_reads or 0),
               min_barcode_reads, min_child_reads))
    #      |----   substep save the raw counts and the reads kept by the clustering for S5 (raw_stats.json)
    os.makedirs(meta_path(args), exist_ok=True)
    with open(meta_path(args, "raw_stats.json"), "w") as fo:
        json.dump({"raw_barcodes": n_distinct, "raw_targets": n_raw_targets,
                   "reads_in": n_reads, "reads_clustered": n_mapped_reads or 0}, fo)

    # step8 : per barcode-cluster statistics (used for the final barcode file)
    #   count      = reads of all barcodes merged into the cluster (coverage after clustering)
    #   n_cluster  = number of different raw barcodes in the cluster
    #   first_cov  = reads of the most abundant barcode, second_cov = of the second one (0 if none)
    # Clusters with a heavy centre are grouped by integer id (narrow rows, no strings; the two
    # largest counts by max(n, 2), no window sort); light barcodes that are their own centre are
    # clusters of size 1.
    #      |----   substep meta/barcode_clusters.parquet: count, n_cluster, first_cov, second_cov per centre
    con.execute("""
        COPY (
            SELECT c.b AS barcode, g."count" AS "count", g.n_cluster AS n_cluster,
                   g.first_cov AS first_cov, g.second_cov AS second_cov
            FROM (SELECT root_id, SUM(n)::BIGINT AS "count", COUNT(*)::BIGINT AS n_cluster,
                         MAX(n)::BIGINT AS first_cov,
                         COALESCE((max(n, 2))[2], 0)::BIGINT AS second_cov
                  FROM (SELECT r.root_id AS root_id, h.n AS n
                        FROM heavy h JOIN hroot r ON r.hid = h.hid
                        UNION ALL
                        SELECT r.root_id AS root_id, m.n AS n
                        FROM lmatch m JOIN hroot r ON r.hid = m.parent_id)
                  GROUP BY root_id) g
            JOIN heavy c ON c.hid = g.root_id
            UNION ALL
            SELECT b, n::BIGINT, 1::BIGINT, n::BIGINT, 0::BIGINT FROM lown
        ) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)
    """ % sql_str(meta_path(args, "barcode_clusters.parquet")))
    log_sub("barcode cluster statistics written")
    #      |----   substep optional barcode -> centre table (--barcode_map)
    if args.barcode_map:
        con.execute("COPY (SELECT barcode, n AS reads, centre AS centre_barcode "
                    "FROM read_parquet(%s)) TO %s (DELIMITER '\t', HEADER)"
                    % (sql_str(bcmap), sql_str(args.out + ".barcode_map.tsv")))
    con.close()

    # step9 : clean up and write the --resume marker
    #      |----   substep delete the temporary DuckDB file and the census Parquet files
    for f in (work + "/bc.duckdb", work + "/bc.duckdb.wal", census_base, census_q):
        if os.path.exists(f):
            os.remove(f)
    #      |----   substep marker file used by --resume
    with open(work + "/bcmap.done", "w") as fo:
        fo.write("ok\n")
