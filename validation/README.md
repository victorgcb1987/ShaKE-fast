# Validation

Ad hoc scripts for proving specific pieces of the pipeline behave as designed,
using synthetic data with a known, hand-computable correct answer - as
opposed to `legacy_scripts/` (old predecessor pipeline entrypoints) or
`src/` (the pipeline itself). See `METHODOLOGY.md` for the biological/design
rationale these scripts are checking against.

## `make_synthetic_hetkmer_dataset.py`

Proves the native het-kmer detection/merging (`src/kmerdata.py`,
`merge_hetkmers` in `src/pipeline.py`) actually works, end to end through the
real `kmc`/`kmc_tools` binaries - not just as isolated unit checks on
hand-built arrays.

It builds one small FASTA file containing four independent, deliberately
constructed loci, each with a known expected outcome:

| Locus | Construction | Expected outcome |
|---|---|---|
| `balanced_snp` | one SNP, allele depths 20 vs 20 | merges in **both** modes (component size 2, coverage ratio ≈1.0) |
| `imbalanced_snp` | one SNP, allele depths 20 vs 1 | merges **only** with `--pooled` (ratio 0.05 fails the non-pooled coverage check) |
| `repeat_hub` | one "hub" k-mer + 7 single-substitution "leaf" neighbors (mimics a repeat/paralog family) | **never** merges (component size 8 exceeds both the pooled cap of 6 and the non-pooled cap of 2) |
| `invariant` | one plain locus, no variants | negative control - no het-kmer edges at all, must be left untouched |

### Usage

```
python validation/make_synthetic_hetkmer_dataset.py -o /tmp/hetkmer_test
python SHaKE.py -i /tmp/hetkmer_test/file_of_files.tsv -o /tmp/hetkmer_test/run_default
python SHaKE.py -i /tmp/hetkmer_test/file_of_files.tsv -o /tmp/hetkmer_test/run_pooled --pooled
```

The generator prints the exact numbers you should see. Check them against:

- **The run log's STEP 4 line** (`grep "components_found" <output_dir>/Kmer_counting_*.log`):
  - default: `components_found=3 components_merged=1 components_rejected=2`
  - `--pooled`: `components_found=3 components_merged=2 components_rejected=1`
- **`Subgroup_Universe_Size` in `results.tsv`**: 12 (default) vs. 11 (`--pooled`), both down from a raw 13 distinct k-mers before any merging.

Last confirmed working (both modes, exact match to predicted numbers, and
confirmed no `.dump`/`.binary`/scratch files left on disk anywhere under the
output directory): 2026-09-21.

If you change the pooled/non-pooled constants in `src/kmerdata.py`
(`POOLED_MAX_COMPONENT_SIZE`, `NONPOOLED_MAX_COMPONENT_SIZE`,
`NONPOOLED_MIN_COVERAGE_RATIO`), rerun this and update the expected numbers
above (and in the script's own printed output) accordingly.
