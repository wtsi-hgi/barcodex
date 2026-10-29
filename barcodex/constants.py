"""
Constants shared by several stages.
"""

# S2 (SQL costs) and S4 (band_dp / cluster_group): one integer carries both edit kinds
COST_SCALE = 1 << 16        # edit cost = hard_edits * COST_SCALE + low_quality_edits
# S1: low 32 bits of an aggregation value = read count
MASK32 = 0xFFFFFFFF
# S1: one 64-bit word of the target low-quality mask (flush_part)
MASK64 = (1 << 64) - 1
# S1: width of the target mask inside the combined mask
TARGET_MASK_BITS = 256      # target low-quality mask width (targets are <= 255 bp)
# S1: longest accepted target
MAX_TARGET_LEN = 255
# S1 + S2: longest barcode that can use quality (2 bits/base in one UBIGINT)
MAX_QUAL_BARCODE_LEN = 32   # barcode quality masks are packed in one 64-bit word (2 bits/base)
# S2: shortest pigeonhole seed
MIN_SEED_LEN = 4            # shortest pigeonhole seed we accept
# S1: chunks buffered per worker queue
QUEUE_DEPTH = 3             # chunks buffered per worker queue
