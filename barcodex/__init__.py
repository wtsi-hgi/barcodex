"""
barcodex: A tool of association of barcodes and targets for paired-end sequencing data.

1. create barcode and target pairs
2. cluster barcodes by sequence similarity and base qualities
3. cluster targets by sequence similarity and base qualities
4. run a statistical test to find significant barcode-target associations
5. fast and memory-efficient implementation using DuckDB and SQL queries
"""

__version__ = "0.1"
