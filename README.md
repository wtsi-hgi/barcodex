# barcodex  (version 0.1)

barcodex links barcodes to the sequences they tag (targets: enhancers, promoters, splicing
variants, ...) in very large paired-end data sets (>10^9 read pairs). It clusters both barcodes
and targets, quality-aware, to remove sequencing and PCR errors. No step keeps all barcodes or
reads in Python memory.

Final map rules: one barcode -> exactly ONE target; one target -> MANY barcodes.

## Install and dependencies

* Python >= 3.10 (tested with 3.10.10 and 3.11.2).
* duckdb >= 1.2.2, < 2 (tested with 1.2.2). Everything else is the Python standard library.
* pigz (optional, tested with 2.6) speeds up reading `.fastq.gz`. Without it, or with
  `--no_pigz`, Python's gzip is used.

The details are in `pyproject.toml`. On the farm, the HGI module has everything:

```bash
module load HGI/softpack/users/fs18/nf_splicing_py
```

`pip install .` in this directory installs the package and a `barcodex` command.

## Running it

```bash
python run_barcodex.py --barcode_read R2.fastq.gz --target_read R1.fastq.gz --out PREFIX \
    --method none none --threads 40 --mem_gb 80

python -m barcodex ...                        # the same, from this directory
python -m barcodex --version
python -m barcodex --selftest                 # synthetic data, all steps, pass / fail
python -m unittest discover -s tests          # unit tests, one file per step
```

All options are long options (`--barcode_read`, `--threads`, ...); there are no short forms.

Every run (and every `--resolve_only`) appends its command to `cmd.log` in the output directory.
Each entry has the date, version, host, working directory, and the command with one option per
line, so it can be copied and run again.

## Choosing how barcode and target are extracted: `--method BARCODE TARGET`

`--method` is required. It takes two values: first the barcode method, then the target method.
They are independent: for example, `--method pos flank` takes the barcode by position and the
target between flanks, and `--method flank none` uses the whole target read.

| method  | the sequence is | options (`*` = `barcode` or `target`) |
|---------|-----------------|---------------------------------------|
| `none`  | the whole read | none |
| `pos`   | a fixed window | `--*_start` (0-based), `--*_len` (0 = to the end of the read) |
| `flank` | between two flanking sequences | `--*_up`, `--*_down`, and `--flank_mismatch` (shared by both, default 2) |

* **pos**: a barcode read that is too short is rejected (`barcode_too_short`). A short target is
  kept as it is.
* **flank**: the flanks are matched on the sequence only, without base qualities, with up to
  `--flank_mismatch` Hamming mismatches. The downstream flank is searched after the upstream
  flank. An exact flank match is used first; otherwise the first window within the mismatch
  limit. This gave the same barcode / target / reject reason for every read as the reference
  scripts: 100,000 reads against `associate_sequences.py` and 200,000 reads against
  `associate_variant_barcode.py`.
* An option of another method is an error, for example `--method none none --barcode_len 16`.

Filters, for any method:

* `--*_min_len N` rejects shorter sequences, for example `--target_min_len 170`.
* `--*_allowed_lens 19,20,21` (or `19-21`) keeps only these lengths.
* `--barcode_temp NNNNATNN...` applies to barcodes only. The barcode must have the template length,
  and its fixed (non-N) bases may differ in at most `--barcode_temp_mismatch` positions
  (default 1). This is the same rule as `check_barcode` in the splicing pipeline.
* Targets longer than 255 bp are rejected (`target_too_long`). This is the width of the target
  quality mask.

With the flank method, use one of these filters for barcodes. Very short barcodes stop S2 with
"barcodes too short".

Clustering uses the base qualities of the extracted sequences. To turn this off, set
`--barcode_max_lowq 0` or `--target_max_lowq 0`.

### Examples

Promoter library (STARR-seq), with the same settings as `run_asso_promoter.sh`:

```bash
python run_barcodex.py --barcode_read Adept_SS_L1_R2.fastq.gz --target_read Adept_SS_L1_R1.fastq.gz \
    --out GP_input --method flank none \
    --barcode_up GTCCTCTTCCCATCT --barcode_down TGGTGGAAGCTACTGTACACCA --barcode_allowed_lens 19,20,21 \
    --threads 40 --mem_gb 80
```

Splicing library combi_6 (merged reads, so both roles read the same file), with the same
settings as `associate_variant_barcode.py` in `run_bc.sh` (see `combi_6_test/run_combi_6.sh`):

```bash
python run_barcodex.py --barcode_read combi_6.extendedFrags.fastq.gz \
    --target_read combi_6.extendedFrags.fastq.gz --out combi_6 \
    --method flank flank --flank_mismatch 2 \
    --target_up tccaatgttccaacctacag --target_down tctttgttctttcagTGAAG --target_min_len 170 \
    --barcode_up GCAGAAAGCACAGAAAGGAA --barcode_down CTACTGATTCGATGCAAGCTT \
    --barcode_temp NNNNATNNNNATNNNNATNNNNATNNNNATNNNNATNN --barcode_temp_mismatch 1 \
    --min_pair_reads 3 --threads 40 --mem_gb 70
```

| `associate_variant_barcode.py` | barcodex |
|---|---|
| `--read1` (variant) / `--read2` (barcode) | `--target_read` / `--barcode_read` |
| flanks for both | `--method flank flank` |
| `--variant_up` / `--variant_down` | `--target_up` / `--target_down` |
| `--variant_check --variant_len 170` | `--target_min_len 170` |
| `--barcode_check --barcode_temp T` | `--barcode_temp T` |
| `--barcode_mismatch 1` | `--barcode_temp_mismatch 1` |
| `--max_mismatches 2` | `--flank_mismatch 2` |
| `--min_barcov 3` | `--min_pair_reads 3` (applied after the clustering) |

## Steps

| step | module | what it does |
|------|--------|--------------|
| S1 ingest   | `fastq.py`, `extract.py`, `ingest.py`, `compaction.py` | One reader process per FASTQ role passes chunks to the parse workers. The workers extract the barcode and the target, count rejects, and aggregate (barcode, target) -> count + low-quality masks in a bounded dict. The dict is written to gz parts, which a background DuckDB compactor merges into Parquet. |
| S2 barcodes | `barcodes.py` | DuckDB barcode census. Barcodes with enough reads ("heavy" barcodes) can be parents. Each barcode is matched to a parent with a pigeonhole seed join and a quality-aware Hamming or Levenshtein distance. Chains of parents are followed to the root by pointer jumping. |
| S3 shards   | `shards.py` | Replaces each barcode by its cluster centre and writes Hive-partitioned shards by hash(centre). |
| S4 targets  | `distance.py`, `targets.py` | Clusters the targets of each barcode, quality-aware, with one process per shard. |
| S5 resolve  | `resolve.py`, `stats.py` | Chooses one target per barcode and writes the result files and the stats. `--resolve_only` re-runs only S5. |

Other modules:

* `cli.py`: options, validation, `cmd.log`, `main`
* `steps.py`: order of the steps and `--resume`
* `common.py`: logging, SQL quoting, resources
* `constants.py`
* `selftest.py`

A mismatch is **hard** when both bases are confident (Phred >= `--*_min_qual`, ACGT); it counts
towards `--*_dist`. It is **low** when one of the two bases has low quality or is N; it counts
towards `--*_max_lowq`. A base is low quality only if it was low quality in every read that
carried that sequence.

## Output (tab-separated, one header line)

`<prefix>.barcode_clusters.tsv` has one row per final barcode:

| column | meaning |
|--------|---------|
| `barcode`, `target` | barcode (centre of its cluster) and its target |
| `barcode_reads` | reads of the barcode cluster (after clustering) |
| `n_diff_barcodes` | number of different raw barcodes merged into the cluster (before clustering) |
| `first_barcode_cov`, `first_barcode_ratio` | reads of the most abundant raw barcode, and its fraction of `barcode_reads` |
| `second_barcode_cov`, `second_barcode_ratio` | the same for the second most abundant raw barcode (0 if none) |

`<prefix>.target_clusters.tsv` has one row per final (barcode, target) pair:

| column | meaning |
|--------|---------|
| `barcode`, `target` | as above |
| `target_reads` | reads of this barcode that support the target (after target clustering) |
| `n_diff_targets` | number of different target sequences merged into it (before clustering) |
| `first_target_cov`, `first_target_ratio` | reads of the most abundant target sequence, and its fraction of `target_reads` |
| `second_target_cov`, `second_target_ratio` | the same for the second most abundant target sequence |

`<prefix>.stats.tsv` has the columns `step`, `stat`, `value`, `percent` and `description`. The
rows are ordered by step:

| step | stats |
|------|-------|
| `1_input` | `total_read_pairs` |
| `2_extraction` | one row per reject reason, then `extracted_read_pairs` |
| `3_barcode_clustering` | `raw_barcodes`, `barcode_clusters`, `reads_dropped_barcode_clustering`, `reads_in_barcode_clusters` |
| `4_target_clustering` | `raw_targets`, `barcode_target_variants` (before) and `barcode_target_clusters` (after) |
| `5_target_assignment` | `barcodes_unique` / `_multi_target` / `_resolved` / `_ambiguous` / `_low_coverage`, and the reads not used: `reads_low_coverage_barcodes`, `reads_ambiguous_barcodes`, `reads_minor_targets` |
| `6_final_map` | final barcodes, targets and pairs, `effective_reads`, mean/median coverage, `mean_barcodes_per_target` |

* `percent` is filled for read counts (`value / total_read_pairs x 100`) and is `NA` otherwise.
* The read counts add up:
  * `total_read_pairs` = the reject rows of `2_extraction` + `extracted_read_pairs`
  * `extracted_read_pairs` = `reads_dropped_barcode_clustering` + `reads_low_coverage_barcodes`
    \+ `reads_ambiguous_barcodes` + `reads_minor_targets` + `effective_reads`
* Reject rows appear only for the reasons that the chosen settings can produce:
  `*_upstream_not_found`, `*_downstream_not_found`, `barcode_too_short`, `*_empty`,
  `target_too_long`, `*_below_min_len`, `*_length_filtered` and `barcode_template_mismatch`.
* Each rejected pair is counted once, under the first reason that applies. The checks run in the
  same order as in the `associate_*.py` scripts: target extraction, barcode extraction, target
  filters, barcode filters.

Also in the output directory:

* `cmd.log`: the commands of all runs.
* `<prefix>.meta/`: the intermediate files for `--resolve_only`. You can delete it, or use
  `--clean_meta`.
* Optional: `<prefix>.barcode_map.tsv` (`--barcode_map`) and `<prefix>_barcode.fastq` /
  `<prefix>_target.fastq` (`--fastq`).

## Tests

| test file | what it covers |
|-----------|----------------|
| `test_extract.py` | none / pos / flank extraction, length and template filters |
| `test_cli.py` | options, `--method`, validation, `cmd.log`, `--version` |
| `test_s1_ingest.py` | FASTQ reader, read parsing and reject order, part files, compaction, step S1 |
| `test_s2_s3_barcodes_shards.py` | barcode clustering (S2) and sharding (S3) |
| `test_s4_targets.py` | edit distance, target clustering within a barcode, step S4 |
| `test_s5_resolve.py` | one target per barcode, output columns, stats file by step |
| `test_all_steps.py` | all steps on synthetic data (`pos none` and `flank pos` scenarios) |

## Code conventions

* Every function has a docstring with a description, Parameters and Returns.
* Steps inside a function are marked `# stepN : text`, and their parts
  `#      |----   substep text`.
* Options are long snake_case options with the same name as the variable:
  `--max_reads` -> `args.max_reads`.
