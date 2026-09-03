# Changelog

All notable changes to this project are documented in this file. See `PERFORMANCE.md` for the earlier performance-focused fixes this fork started from.

## [Unreleased]

### Changed
- Removed the `--presence`/`-p` CLI flag. Instead of choosing between count-based and presence/absence diversity calculation, `SHaKE.py` now always computes both and reports them side by side in `results.tsv`, which gained five new columns: `Diversity_log2_presence`, `Specifity_log2_presence`, `Diversity_log10_presence`, `Specifity_log10_presence`, `Kolmogorov_presence`. Applies to genomic/transcriptome and expression samples alike.
- Refactored `SHaKE.py`: pipeline orchestration (previously a ~270-line sequential `main()`) moved into a new `src/pipeline.py`, split into named stage functions (`build_databases`, `build_histograms`, `dump_counts`, `merge_hetkmers`, `compute_universe_sizes`, `compute_estimators`, `compute_expression_estimators`, `write_outputs`, `run_pipeline`). `SHaKE.py` is now a thin CLI entrypoint. The existing performance-optimized bodies of `src/kmer.py`, `src/kolmogorov.py`, and `src/utils.py` (single-pass entropy calc, batched I/O, merge-sort universe size — see `PERFORMANCE.md`) are unchanged; only their signatures/return values gained the plumbing described below.
- Extracted the inlined het-kmer union-find merge logic into `src/kmc.py:merge_kmers_by_hetkmers`, following the `check_run()` convention.
- `src/kmer.py:calculate_sample_shannon_estimators` gained a `suffix` parameter and now merges into an existing `estimators[group][sub][name]` entry instead of overwriting it, so a regular and a presence/absence pass over the same sample can land in the same row (needed for the presence/absence change above).
- `src/kolmogorov.py:calculate_kolmogorov_estimator`, `create_kmer_binary_file`, and `create_expression_binary_file` now return a status dict (command/returncode/msg/out_fpath) instead of `None`/a bare path, so callers can log their outcome; `calculate_kolmogorov_estimator` also gained a `key` parameter (default `"kolmogorov"`) so a regular and a presence pass can both be recorded without clobbering each other.
- `src/expression.py:calculate_sample_estimators` now also computes `diversity_log2`/`specifity_log2` (previously log10-only); existing `diversity`/`specifity` keys are unchanged for backward compatibility with `legacy_scripts/calculate_expression_diversity_by_sample.py`.
- Added `src/utils.py:log_and_print` to collapse the repeated print/write/flush logging pattern.

### Fixed
- `results.tsv` rows for `expression`-kind samples reported the diversity/specificity/kolmogorov values of the last genomic sample processed instead of their own, and corrupted an unrelated sample's stored kolmogorov value. Expression rows now compute and report their own values correctly.
- `create_expression_binary_file` never actually received or used the `presence` flag — expression-kind Kolmogorov values were always computed from raw TPM values regardless of presence/absence mode, unlike genomic/transcriptome data where presence encoding worked correctly. Now `presence` is passed through and honored for expression data too.
- `file_manifiest.tsv` header was malformed (`Group\tSubgroup\tRep\tKind"File`, missing a tab, stray quote). Fixed to `Group\tSubgroup\tRep\tKind\tFile`, matching README.md.
- `main()` in `SHaKE.py` called `get_arguments()` twice, discarding the first call's result (leaking an unclosed log file handle) and silently re-running STEP 0 (BAM→FASTA conversion) twice. Now called once.
- The het-kmer union-find merge step and `create_kmer_binary_file` skipped already-completed work silently, with no status line. Both now report status through `check_run()`, so reruns log `#ALREADY_DONE:` consistently across all stages.

### Notes
- No change to the CLI flags (besides removing `--presence`), input file-of-files format, or intermediate file naming/paths — an output directory from a prior run remains resumable.
- Known, intentionally unfixed: `dump_counts` reads `lowerbound`/`upperbound` from the wrong dict level, so the dump step always falls back to its defaults (1 / ~10 billion) regardless of the input file's bound columns. STEP 1's kmc counting already enforces the real bounds via `-ci`/`-cx`, so this is likely harmless double-filtering, but fixing it would change dump-file contents — left for a separate, explicitly-scoped change.
