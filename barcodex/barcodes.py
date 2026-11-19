"""
Step 2 -- processing barcodes: 
Clustering barcodes by sequence similarity and abundance, then mapping every barcode to its cluster centre.

1. The parental barcodes are the ones with count >= --barcode_ratio (the "heavy" set). 
2. Every barcode is matched to its parents by pigeonhole-seed JOIN + quality-aware Hamming (or Levenshtein for indels).

SQL builders: packed barcodes, seed join, ranking, matching queries
"""

import glob
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
            bit_shift = f"{fixed_len - 1 - i}"
            lane_shift = f"{2 * (fixed_len - 1 - i)}"
        else:
            bit_shift = f"greatest(length(b) - 1 - {i}, 0)"
            lane_shift = f"(2 * {bit_shift})"
        #    |---- substep code lane: A=0 C=1 G=2 T=3 (anything else = 0)
        code = (
            f"(CASE substr(b, {i + 1}, 1) "
            "WHEN 'C' THEN 1 "
            "WHEN 'G' THEN 2 "
            "WHEN 'T' THEN 3 "
            "ELSE 0 END)"
            f"::UBIGINT << ({lane_shift})::UBIGINT"
        )
        #    |---- substep confidence lane: 1 when the mask bit of the base is 0 (confident)
        conf = (
            "(1::UBIGINT - "
            f"((bm >> ({bit_shift})::UBIGINT) & 1::UBIGINT))"
            f" << ({lane_shift})::UBIGINT"
        )
        #    |---- substep variable length: positions beyond the barcode end contribute 0
        if not fixed_len:
            code = (
                f"CASE WHEN length(b) > {i} "
                f"THEN {code} "
                "ELSE 0::UBIGINT END"
            )
            conf = (
                f"CASE WHEN length(b) > {i} "
                f"THEN {conf} "
                "ELSE 0::UBIGINT END"
            )
        code_terms.append(f"({code})")
        conf_terms.append(f"({conf})")

    # step2 : OR all lane terms together into two 64-bit words
    return (
        f"({' | '.join(code_terms)})::UBIGINT",
        f"({' | '.join(conf_terms)})::UBIGINT"
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
        from_clause = (
            f"{src} c CROSS JOIN range({n_seeds}) AS s(slot) "
            f"CROSS JOIN range({-band}, {band + 1}) AS h(sh)"
        )
        pos = f"s.slot * {seed_len} + h.sh"
        shift_filter = (
            f"AND {pos} >= 0 "
            f"AND {pos} + {seed_len} <= length(c.b)"
        )
        seed_expr = f"hash(substr(c.b, {pos} + 1, {seed_len}))"
        length_filter = f"AND abs(length(c.b) - length(p.b)) <= {band}"
    # without indels: only the exact slot positions, no shift / length filters
    else:
        from_clause = f"{src} c CROSS JOIN range({n_seeds}) AS s(slot)"
        seed_expr = f"hash(substr(c.b, s.slot * {seed_len} + 1, {seed_len}))"
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
        f"({cost_expr})::HUGEINT * 79228162514264337593543950336 + "                                
        f"(4294967295 - LEAST({parent_count_expr}, 4294967295))::HUGEINT * 18446744073709551616 + " 
        f"({tiebreak_expr})::HUGEINT"                                                               
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
        dist = (
            "CASE WHEN length(c.b) = length(p.b) "
            "THEN hamming(c.b, p.b) "
            "ELSE 99 END"
        )

    # ranking key for arg_min: smaller distance, then larger parent count, then hash
    key = rank_key_sql(dist, "p.n", "hash(p.b)")

    # assemble: seed join -> heavy lookup -> coverage / tie rules -> distance filter -> best parent
    return (
        "SELECT c.b AS barcode, "
        "c.n AS n, "
        f"arg_min(p.hid, {key}) AS parent_id "
        f"FROM {from_clause} "
        f"JOIN hseed q ON q.slot = s.slot AND q.seed_hash = {seed} "
        "JOIN heavy p ON p.hid = q.hid "
        "WHERE p.b <> c.b "
        f"AND p.n >= {float(ratio)!r} * c.n "
        "AND (p.n > c.n OR (p.n = c.n AND p.b < c.b)) "
        f"{shift_filter} "
        f"{length_filter} "
        f"AND ({dist}) <= {max_dist} "
        "GROUP BY c.b, c.n"
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
        "bit_count(("
        f"(({xor_expr}) | (({xor_expr}) >> 1::UBIGINT))"    # mismatching bases (either lane bit differs)
        " & c.conf & p.conf"                                # ... where both bases are confident
        ")::BIGINT)"
    )

    lev_col = ", levenshtein(c.b, p.b) AS lev" if band > 0 else ""

    # Step 3: Build candidate pairs and calculate their mismatch statistics.
    # Each candidate must:
    #   - Have a different barcode sequence from the child.
    #   - Meet the minimum parent-to-child read-count ratio.
    #   - Follow the abundance-based ordering used to avoid duplicate child-parent comparisons.
    layer_a = (
        "SELECT c.b AS barcode, "
        "c.n AS n, "
        "p.hid AS parent_id, "
        "p.n AS parent_n, "
        "hash(p.b) AS parent_hash, "
        f"CASE WHEN {same_len} THEN hamming(c.b, p.b) END AS total, "
        f"CASE WHEN {same_len} THEN {hard_expr} END AS hard"
        f"{lev_col} "
        f"FROM {from_clause} "
        f"JOIN hseed q ON q.slot = s.slot AND q.seed_hash = {seed} "
        "JOIN heavy p ON p.hid = q.hid "
        "WHERE p.b <> c.b "
        f"AND p.n >= {float(ratio)!r} * c.n "
        "AND (p.n > c.n OR (p.n = c.n AND p.b < c.b)) "
        f"{shift_filter} "
        f"{length_filter}"
    )

    # Step 4: Evaluate candidate pairs using two independent matching criteria.
    # Hamming cost:
    #   Accept only if both the high-quality mismatch count and the low-quality mismatch count are within limits.
    #   low-quality mismatches = total mismatches - high-quality mismatches
    # Levenshtein cost:
    #   Accept pairs whose edit distance is within max_hard.
    ham_cost = (
        "CASE WHEN total IS NOT NULL "
        f"AND hard <= {max_hard} "
        f"AND total - hard <= {max_low} "
        f"THEN hard * {COST_SCALE} + (total - hard) END"
    )

    if band > 0:
        lev_cost = (
            f"CASE WHEN lev <= {max_hard} "
            f"THEN lev * {COST_SCALE} END"
        )
    else:
        lev_cost = "NULL"

    layer_b = (
        "SELECT barcode, n, parent_id, parent_n, parent_hash, "
        f"{ham_cost} AS ham_cost, "
        f"{lev_cost} AS lev_cost "
        f"FROM ({layer_a})"
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
        "SELECT barcode, "
        "n, "
        f"arg_min(parent_id, {key}) AS parent_id "
        f"FROM ({layer_b}) "
        "WHERE ham_cost IS NOT NULL OR lev_cost IS NOT NULL "
        "GROUP BY barcode, n"
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
        action = "enable" if enabled else "reset"
        log_sub(f"NOTE: could not {action} join-order pinning ({str(e).splitlines()[0]})")


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
        log_sub(f"query plan written to {path}")
    except Exception as e:
        log_sub(f"NOTE: could not write the query plan ({str(e).splitlines()[0]})")


# ---- Stage processing -----------------------------------------------------

def stage_barcodes(args, work, threads, mem_gb):
    """
    S2: Build barcode clusters by assigning each raw barcode to a centre barcode.

    1. Count reads for each distinct barcode and collect its quality information.
    2. Identify heavy barcodes, which are eligible to act as cluster parents.
    3. Match each barcode to its best heavy parent using seed-based matching.
       Matching uses quality-aware Hamming distance, or Levenshtein distance when indels are enabled.
    4. Resolve parent chains among heavy barcodes so that every barcode points directly to its final cluster centre.
    5. Assign light barcodes to the centre of their matched heavy parent.
       Unmatched light barcodes are retained as their own centres only when they meet the minimum read-count threshold.
    6. Calculate barcode-cluster statistics and write the output files.

    Parameters:
        -- args: argparse.Namespace
        -- work: working directory for intermediate and output files
        -- threads: DuckDB threads
        -- mem_gb: RAM budget in GB (85% is given to DuckDB)
    Returns:
        -- None: results are files
    """
    # step1: initialise working files and DuckDB
    census_base = work + "/census_base.parquet"
    census_q = work + "/census_q.parquet"
    bcmap = work + "/bcmap.parquet"
    tuples = f"read_parquet({sql_str(work + '/tuples/t*.parquet')})"
    log("S2: barcode census + clustering (DuckDB, out-of-core) ...")
    log_sub(f"DuckDB budget for S2: {mem_gb * 0.85:.1f} GB RAM, {threads} thread(s); "
            f"spill directory {work}/duck_bc")
    con = open_duckdb(work, threads, mem_gb * 0.85, "bc", work + "/bc.duckdb")

    # step2: Build the barcode census.
    con.execute(
        "COPY (SELECT barcode AS b, SUM(cnt)::BIGINT AS n, bit_and(bm) AS bm "
        f"FROM {tuples} "
        "GROUP BY barcode) "
        f"TO {sql_str(census_base)} (FORMAT PARQUET, COMPRESSION ZSTD)"
    )

    n_distinct, n_reads, min_len, max_len = con.execute(
        "SELECT COUNT(*), SUM(n), MIN(length(b)), MAX(length(b)) "
        f"FROM read_parquet({sql_str(census_base)})"
    ).fetchone()
    log_sub(f"{n_reads} reads, {n_distinct} distinct barcodes (length {min_len}-{max_len})")

    n_raw_targets = con.execute(
        "SELECT COUNT(*) "
        f"FROM (SELECT target FROM {tuples} GROUP BY target)"
    ).fetchone()[0]
    log_sub(f"{n_raw_targets} distinct raw target sequences")

    tuple_files = glob.glob(work + "/tuples/t*.parquet")
    if len(tuple_files) == 1:
        n_raw_pairs = con.execute(f"SELECT COUNT(*) FROM {tuples}").fetchone()[0]
    else:
        n_raw_pairs = con.execute(
            "SELECT COUNT(*) "
            "FROM (SELECT md5_number(barcode || '|' || target) AS k "
            f"FROM {tuples} "
            "GROUP BY k)"
        ).fetchone()[0]
    log_sub(f"{n_raw_pairs} distinct raw (barcode, target) pairs")

    # step3: ratio, distances, quality mode and seed layout
    ratio = max(args.barcode_ratio, 1.0)   
    max_hard = args.barcode_dist      
    band = min(args.barcode_indel, max_hard)
    use_quality = args.barcode_max_lowq > 0 and max_len <= MAX_QUAL_BARCODE_LEN
    max_low = args.barcode_max_lowq if use_quality else 0

    if args.barcode_max_lowq > 0 and not use_quality:
        log_sub(f"NOTE: barcodes longer than {MAX_QUAL_BARCODE_LEN} bp: "
                "quality-aware clustering off for barcodes")

    #      |----   substep pigeonhole seeds: pairs with <= (max_hard + max_low) differences share one of (max_hard + max_low + 1)
    #      |----   substep seeds exactly; lower max_low while the seeds would get shorter than MIN_SEED_LEN
    while use_quality and max_low > 0 and min_len // (max_hard + max_low + 1) < MIN_SEED_LEN:
        max_low -= 1
        log_sub(f"NOTE: barcodes too short for the requested --barcode_max_lowq; using {max_low}")

    n_seeds = max_hard + max_low + 1
    seed_len = min_len // n_seeds
    if seed_len < MIN_SEED_LEN:
        raise SystemExit(f"ERROR: barcodes too short for --barcode_dist {max_hard}")
    
    #      |----   substep read-count thresholds for dropping light barcodes
    min_barcode_reads = (args.min_barcode_reads if args.min_barcode_reads is not None else args.min_pair_reads)
    min_child_reads = max(1, args.min_child_reads)

    # step4: census table used by all later queries
    if use_quality:
        code_expr, conf_expr = packed_barcode_sql(max_len, max_len if min_len == max_len else 0)
        con.execute(
            "COPY (SELECT b, "
            "n, "
            f"{code_expr} AS code, "
            f"{conf_expr} AS conf "
            f"FROM read_parquet({sql_str(census_base)})) "
            f"TO {sql_str(census_q)} (FORMAT PARQUET, COMPRESSION ZSTD)"
        )
        census_file = census_q
        columns = "b, n, code, conf"
        log_sub(f"quality-aware matching: <= {max_hard} hard + <= {max_low} low-quality mismatches, "
                f"{n_seeds} seeds of {seed_len} bp")
    else:
        census_file = census_base
        columns = "b, n"

    def matched(src):
        """
        SQL that matches every row to its best heavy parent
        Parameters:
            -- src: str: SQL table/subquery with columns b, n (+ code, conf when quality-aware)
        Returns:
            -- str: SQL query returning (barcode, n, parent_id)
        """
        if use_quality:
            return matched_query_quality(src, seed_len, n_seeds, ratio, max_hard, max_low, band)
        return matched_query_plain(src, seed_len, n_seeds, ratio, max_hard, band)

    #      |----   substep census source for the following queries
    census = f"read_parquet({sql_str(census_file)})"

    # step5: heavy barcodes = the only possible parents; build their seed index
    # Every heavy barcode gets an integer id (hid).  All big joins below use ids and 64-bit
    # seed hashes instead of the barcode strings, so the join tables stay small and in memory.
    #      |----   substep heavy table: integer id (hid) + packed barcode columns
    con.execute(
        "CREATE OR REPLACE TABLE heavy AS "
        f"SELECT ROW_NUMBER() OVER () AS hid, {columns} "
        f"FROM {census} "
        f"WHERE n >= {float(ratio)!r}"
    )
    n_heavy = con.execute("SELECT COUNT(*) FROM heavy").fetchone()[0]
    log_sub(f"{n_heavy} heavy barcodes (count >= {ratio:g}) can be cluster centres/parents")

    #      |----   substep seed index of the heavy set: (slot, seed hash, heavy id); seeds shared by more than --max_cand
    #      |----   substep heavy barcodes (low complexity) are dropped to bound the join size
    con.execute(f"""
        CREATE OR REPLACE TABLE hseed AS
        SELECT hid, slot, seed_hash FROM (
            SELECT h.hid AS hid, s.slot AS slot,
                   hash(substr(h.b, s.slot * {seed_len} + 1, {seed_len})) AS seed_hash
            FROM heavy h CROSS JOIN range({n_seeds}) AS s(slot))
        QUALIFY count(*) OVER (PARTITION BY slot, seed_hash) <= {args.max_cand}
    """)
    n_seed_rows = con.execute("SELECT COUNT(*) FROM hseed").fetchone()[0]
    log_sub(f"seed index built: {n_seed_rows} entries "
            f"({n_seeds} seeds of {seed_len} bp per heavy barcode); "
            f"spill so far {spill_gb(work):.1f} GB")

    # step6: heavy barcodes: best heavy parent, then pointer jumping to the root
    #      |----   substep pin the join order: children JOIN integer seed index JOIN heavy strings
    pin_join_order(con, True)
    #      |----   substep best heavy parent of every heavy barcode (hroot: hid -> root_id)
    con.execute(
        "CREATE OR REPLACE TABLE hroot AS "
        "SELECT h.hid AS hid, COALESCE(m.parent_id, h.hid) AS root_id "
        "FROM heavy h "
        f"LEFT JOIN ({matched('heavy')}) m ON m.barcode = h.b"
    )
    n_with_parent = con.execute("SELECT COUNT(*) FROM hroot WHERE root_id <> hid").fetchone()[0]
    log_sub(f"{n_with_parent} of {n_heavy} heavy barcodes have a more abundant heavy parent; "
            f"spill so far {spill_gb(work):.1f} GB")
    
    #      |----   substep pointer jumping: root := root of root, until nothing changes
    jump_round = 0
    while True:                    
        jump_round += 1
        con.execute(
            "CREATE OR REPLACE TABLE hroot2 AS "
            "SELECT a.hid AS hid, COALESCE(b.root_id, a.root_id) AS root_id "
            "FROM hroot a "
            "LEFT JOIN hroot b ON a.root_id = b.hid"
        )
        n_changed = con.execute(
            "SELECT COUNT(*) "
            "FROM hroot a "
            "JOIN hroot2 b USING (hid) "
            "WHERE a.root_id <> b.root_id"
        ).fetchone()[0]
        con.execute("DROP TABLE hroot")
        con.execute("ALTER TABLE hroot2 RENAME TO hroot")
        log_sub(f"pointer jumping round {jump_round}: {n_changed} heavy barcodes changed root")
        if not n_changed:
            break

    #      |----   substep centre barcode (string) of every heavy barcode (hcentre)
    con.execute(
        "CREATE OR REPLACE TABLE hcentre AS "
        "SELECT r.hid AS hid, r.root_id AS root_id, c.b AS centre "
        "FROM hroot r "
        "JOIN heavy c ON c.hid = r.root_id"
    )

    # step7: every barcode -> centre barcode
    #  heavy barcode            -> its root
    #  light barcode, parent    -> root of its heavy parent
    #  light barcode, no parent -> its own centre, only if it has >= min_barcode_reads reads
    children = (
        f"(SELECT {columns} "
        f"FROM {census} "
        f"WHERE n < {float(ratio)!r} "
        f"AND n >= {min_child_reads})"
    )

    #      |----   substep save the query plan to check the join order
    dump_plan(con, matched(children), work + "/explain_light_match.txt")
    log_sub("matching light barcodes to heavy parents (the longest step) ...")

    #      |----   substep match light barcodes to heavy parents
    con.execute(
        "CREATE OR REPLACE TABLE lmatch AS "
        "SELECT barcode, n, parent_id, hash(barcode) AS hk "
        f"FROM ({matched(children)})"
    )

    #      |----   substep release the join-order pin
    pin_join_order(con, False)
    n_lmatch = con.execute("SELECT COUNT(*) FROM lmatch").fetchone()[0]
    log_sub(f"{n_lmatch} light barcodes matched to a heavy parent; "
            f"spill so far {spill_gb(work):.1f} GB")

    #      |----   substep light barcodes without a parent keep their own centre when they have enough reads
    #      |----   substep (anti-join on the 64-bit barcode hash: no strings in the join table)
    con.execute(f"""
        CREATE OR REPLACE TABLE lown AS
        SELECT s.b AS b, s.n AS n FROM {census} s
        WHERE s.n >= {max(min_barcode_reads, min_child_reads)} AND s.n < {float(ratio)!r}
          AND NOT EXISTS (SELECT 1 FROM lmatch m WHERE m.hk = hash(s.b))
    """)
    n_lown = con.execute("SELECT COUNT(*) FROM lown").fetchone()[0]
    log_sub(f"{n_lown} light barcodes without a parent kept as their own cluster")

    #      |----   substep bcmap.parquet = heavy -> root | light with parent -> root of parent | light own centre
    con.execute(f"""
        COPY (
            SELECT h.b AS barcode, c.centre AS centre, h.n AS n
            FROM heavy h JOIN hcentre c ON c.hid = h.hid
            UNION ALL
            SELECT m.barcode, c.centre, m.n                             
            FROM lmatch m JOIN hcentre c ON c.hid = m.parent_id
            UNION ALL
            SELECT b, b, n FROM lown                                     
        ) TO {sql_str(bcmap)} (FORMAT PARQUET, COMPRESSION ZSTD)
    """)

    #      |----   substep report mapped / dropped barcodes
    n_mapped, n_centres, n_mapped_reads = con.execute(
        "SELECT COUNT(*), COUNT(*) FILTER (WHERE barcode = centre), SUM(n) "
        f"FROM read_parquet({sql_str(bcmap)})"
    ).fetchone()
    n_dropped_reads = n_reads - (n_mapped_reads or 0)
    log_sub(f"{n_mapped} distinct barcodes -> {n_centres} barcode clusters; "
            f"dropped {n_distinct - n_mapped} barcodes ({n_dropped_reads} reads) "
            f"(unmatched with < {min_barcode_reads} reads, "
            f"or < {min_child_reads} reads = --min_child_reads)")
    
    #      |----   substep save the raw counts and the reads kept by the clustering for S5 (raw_stats.json)
    os.makedirs(meta_path(args), exist_ok=True)
    with open(meta_path(args, "raw_stats.json"), "w") as fo:
        json.dump({"raw_barcodes": n_distinct, "raw_targets": n_raw_targets,
                   "raw_barcode_target_pairs": n_raw_pairs,
                   "reads_in": n_reads, "reads_clustered": n_mapped_reads or 0}, fo)

    # step8: per barcode-cluster statistics (used for the final barcode file)
    #   count      = reads of all barcodes merged into the cluster (coverage after clustering)
    #   n_cluster  = number of different raw barcodes in the cluster
    #   first_cov  = reads of the most abundant barcode, second_cov = of the second one (0 if none)
    clusters_out = sql_str(meta_path(args, "barcode_clusters.parquet"))
    con.execute(f"""
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
        ) TO {clusters_out} (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    log_sub("barcode cluster statistics written")

    #      |----   substep optional barcode -> centre table (--barcode_map)
    if args.barcode_map:
        map_out = sql_str(args.out + ".barcode_map.tsv")
        con.execute(
            "COPY (SELECT barcode, n AS reads, centre AS centre_barcode "
            f"FROM read_parquet({sql_str(bcmap)})) "
            f"TO {map_out} (DELIMITER '\t', HEADER)"
        )
    con.close()

    # step9: clean up and write the --resume marker
    for f in (work + "/bc.duckdb", work + "/bc.duckdb.wal", census_base, census_q):
        if os.path.exists(f):
            os.remove(f)

    with open(work + "/bcmap.done", "w") as fo:
        fo.write("ok\n")
