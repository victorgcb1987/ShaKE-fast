# Shake - Shannon diversity on K-mers and gene Expression

## Requirements

kmc https://github.com/refresh-bio/KMC. You will need to add kmc binaries to PATH variable for example using `export PATH=$PATH:kmc/install/bin/`

**numpy**: used for k-mer encoding and het-kmer detection. `pip install numpy`.

(smudgeplot is no longer required — het-kmer detection is now done natively, see `--pooled` below.)

clone KATULU respository: `git clone https://github.com/victorgcb1987/ShaKE.git`

## How to use 

### Omics diversity calculator
Omics diversity calculator is a pipeline that will calculate Shannon diversity index for files of your choice. These file can be fasta/fastq files or gene expression files like the ones provided by stringtie.

In order to runt this program you will need to prepare a file with table similar to this one:

| onegroup | R | transcriptome | 0 | 9999999999 | test1_r1.fastq.gz,test1_r2.fastq.gz |
|---|---|---|---|---|---|
| onegroup | R | transcriptome | 0 | 9999999999 | test2_r1.fastq.gz,test2_r2.fastq.gz |
| onegroup | G | genome | 0 | 9999999999 | test5.fasta.gz |
| othergroup | G | transcriptome | 0 | 999999999 | test4_r1.fastq.gz,test4_r2.fastq.gz |
| othergroup | E | expression | 0 | 000000000 | Petal.guided.abund.tsv |
| othergroup |E | expression | 0 | 000000000  |adult_vascular_leaf.guided.abund.tsv |
| anothergroup | G  |genome | 0 | 9999999999 | test5.fasta.gz |

The first and second columns are used to label group and subgroup respectively, for example to mark that one tissue in the first column and then a set of files of the same origin/class/etc. with the second column. 

The third column serves to label what kind of data is; you can put whaterver you like but if you use *** transcriptome *** hetkmers will be calculated and kmers and their values will be grouped using this (see `--pooled` below for how this detection behaves). If you put *** expression ***, files will be treated as expression tables.

Fourth and fifth column are the minimun and the maximum values cutoff in order to consider a kmer or not in the analysis.

Sixth column is the filepaths column, you can separate different files with a comma.

This pipeline can be run with `omics_diversity_pipeline.py`and it accepts the following arguments:

  `--input_file INPUT_FILE, -i INPUT_FILE` (Required): the path to the table previously described
  
  `--output_dir OUTPUT_DIR, -o OUTPUT_DIR`(Required): the path were interemediate files and results will be stored. Building each sample's kmc database (the slow step) is skipped on a rerun if it already exists in this directory; the k-mer dump, het-kmer detection, and diversity/kolmogorov estimators are cheap enough now that they're always recomputed from the kmc database rather than cached to disk (see "Output" below - no `.dump`/binary scratch files are written at all anymore).
  
  `--ram_usage RAM_USAGE, -r RAM_USAGE` (Optional, 1 by default): The max memory usage of kmc.
  
  `--num_threads NUM_THREADS, -t NUM_THREADS`(Optional 1 by default): Number of threads used by kmc.
  
  `--kmer_size KMER_SIZE, -k KMER_SIZE` (Optional, 21 by default): the kmer size used by kmc.
  
  `--merge_universe, -m` (Optional, False by default): This switchs whether if the kmer universe size should be calculated from all items from a group (True) or only for each sub group.

  `--pooled` (Optional, False by default): controls how het-kmer detection decides whether two k-mers one substitution apart represent the same locus (e.g. a SNP) rather than an unrelated/repeat-family collision. With `--pooled` (multiple individuals pooled into one sample, so allele frequency at a site can be anything and more than 2 alleles can legitimately co-occur), it only rejects a group of connected k-mers once it gets implausibly large (a repeat-family signature). Without it (the default - a single individual, at most 2 alleles per site), it additionally requires the pair's coverage to be reasonably balanced, closer to what a true heterozygous site looks like in one diploid genome. Het-kmer detection is done natively (no external tool) either way.

  `--sequential, -s` (Optional, False by default): low-disk/low-memory mode. kmc's `.kmc_pre`/`.kmc_suf` databases are what take the space, and normally all of them (and all the k-mer tables) are kept until the end of the run. With this flag each sample goes through its whole chain on its own (kmc database, histogram, in-memory k-mer table, het-kmer merge) and is then reduced to what the results need, and its database and table are dropped before the next sample starts. The disk holds one database at a time and the memory one k-mer table at a time. `results.tsv` and `file_manifiest.tsv` are written at the end and are identical to a normal run. What a sample keeps meanwhile: its Shannon diversity (it doesn't depend on the universe size) and a streaming Kolmogorov compressor already fed with its own counts. Specificity and the Kolmogorov ratio depend on the universe size, so they are finished when the last sample of the group is done; the only extra thing held until then is, with `--merge_universe`, one running union of the k-mers defining the group's universe. Notes: only databases built by that run are deleted (a pre-existing one is left alone); the `.hist` and `.files` files are kept; since the databases are gone, rerunning recounts every sample.

  Every run also computes a presence/absence variant of the diversity/kolmogorov estimators alongside the regular count-based ones, reporting both in the same `results.tsv` row (see below) — counting all kmer count values to 1 if present (or, for expression tables, a gene gets a value of 1 if its TPM is 1 or more, 0 otherwise). There is no separate flag for this; it's always calculated.

  ## Output

  ### Diversity results

  Shake will generate a table called `results.tsv` inside the directory specified with the `--output_dir` argument. This table has the following format:

| Group      | Subgroup | Rep              | Kind         | Subgroup_Universe_Size | Diversity_log2 | Specifity_log2 | Diversity_log10 | Specifity_log10 | Kolmogorov         | Diversity_log2_presence | Specifity_log2_presence | Diversity_log10_presence | Specifity_log10_presence | Kolmogorov_presence |
|------------|----------|------------------|--------------|------------------------|----------------|----------------|------------------|------------------|---------------------|--------------------------|---------------------------|----------------------------|-----------------------------|-----------------------|
| ERR3333388 | R        | ERR3333388_R1    | transcriptome| 40072639               | 22.7620354967  | 2.4940786897   | 6.8520554469     | 0.7507924971     | 0.0636632271        | *(same 5 metrics, presence/absence)* | | | | |
| ERR3333389 | R        | ERR3333389_R1    | transcriptome| 38691164               | 22.9552288052  | 2.2502719911   | 6.9102124277     | 0.6773993677     | 0.0648244898        | *(same 5 metrics, presence/absence)* | | | | |
| ERR3333390 | R        | ERR3333390_R1    | transcriptome| 39507337               | 22.9614172079  | 2.2742000608   | 6.9120753225     | 0.6846024344     | 0.0653194914        | *(same 5 metrics, presence/absence)* | | | | |
| ERR3333391 | R        | ERR3333391_R1    | transcriptome| 40188471               | 22.5576212642  | 2.7026570900   | 6.7905206314     | 0.8135808521     | 0.0647497717        | *(same 5 metrics, presence/absence)* | | | | |
| ERR3333392 | R        | ERR3333392_R1    | transcriptome| 36116326               | 21.8123101535  | 3.2938376502   | 6.5661596309     | 0.9915439335     | 0.0630993513        | *(same 5 metrics, presence/absence)* | | | | |

(the columns after Kolmogorov are the real, currently-generated `Diversity_log2_presence`/`Specifity_log2_presence`/`Diversity_log10_presence`/`Specifity_log10_presence`/`Kolmogorov_presence` columns; sample values are omitted here since this table predates that addition)

  `Group , Subgroup, Rep and Kind`: have the values provided by the user in the file of files.
  
  `Subgroup_Universe_Size`: is the number of unique K-mers (universe size) found in the dataset.
  
  `Diversity_log2 and Diversity_log10`: is the Shannon Diversity Index score calculated to for each dataset, using log2 and log10, from the raw K-mer/expression counts.
  
  `Specificty_log2 and Specifity_log10`: is the K-mer specifity of each dataset. A specifity value of 0 means that each kind of K-mer is equally represented in a given dataset and higher values means that some K-mers are overrepresented. Specificty is calculated for log2 and log10.
  
  `Kolmogorov`: is Kolmogorov complexity measurement. Is an alternative measurement for Shannon Diversity Index that reduces the impact of K-mers exclusive to unique datasets from the same group.

  `Pielou_Evenness` and `Universe_Evenness`: Shannon diversity rescaled to 0-1 (the ratio is the same in any log base). `Pielou_Evenness = H / log(n)`, with `n` the number of k-mers (genes, for expression) actually present in the sample: it only depends on the sample itself, so it doesn't change when other samples are added or removed. `Universe_Evenness = H / log(Subgroup_Universe_Size)` (equivalently `1 - Specifity / log(universe)`): it also reflects how much of the universe the sample covers, so it depends on the other samples of the subgroup/group (or group with `--merge_universe`). Both are `nan` when the denominator is 0 (a single k-mer or universe of 1) and are not clipped, so `Universe_Evenness` can marginally exceed 1 if a sample has more k-mers than its universe (possible after het-kmer merging). Computed from the regular, count-based diversity only. For expression samples the universe is the number of genes in the table.

  `Kolmogorov_norm`: the Kolmogorov ratio placed between two references computed with the same number of k-mers and the same universe size: `(C - floor) / (reference - floor)`, where `C` is the compressed size of the sample's count profile, `floor` the compressed size of the presence/absence profile (the least complex profile for that n and universe) and `reference` the compressed size of i.i.d. geometric counts with the sample's own mean count (the maximum-entropy distribution over non-negative integers for that mean, generated with a fixed seed so it is reproducible). Roughly 0 means a profile as regular as presence/absence, about 1 means as unpredictable as a maximum-entropy profile of the same depth; it is not clipped, so it can come out slightly above 1, and it is `nan` when undefined (no k-mers). Unlike the raw `Kolmogorov` it does not depend on the sequencing depth or on how large the universe is. It is `NA` for expression samples and there is no presence/absence version (that profile only depends on n and the universe size).

  `Diversity_log2_presence, Specifity_log2_presence, Diversity_log10_presence, Specifity_log10_presence and Kolmogorov_presence`: the same five metrics recalculated using presence/absence instead of raw counts (every present K-mer/gene counts as 1 regardless of its actual count/TPM). Specifity values here will be 0 or near 0, since presence/absence flattens out over/under-representation. These are always computed alongside the regular columns — there is no separate flag to enable them.

  ### File Manifiest
  Other file found in output dir will be `file_manifiest.tsv`with the following format:

  | Group      | Subgroup | Rep            | Kind         | File                              |
|------------|----------|----------------|--------------|-----------------------------------|
| ERR3333388 | R        | ERR3333388_R1  | transcriptome| ERR3333388_q30l50_R1.fastq.gz     |
| ERR3333388 | R        | ERR3333388_R1  | transcriptome| ERR3333388_q30l50_R2.fastq.gz     |
| ERR3333389 | R        | ERR3333389_R1  | transcriptome| ERR3333389_q30l50_R1.fastq.gz     |
| ERR3333389 | R        | ERR3333389_R1  | transcriptome| ERR3333389_q30l50_R2.fastq.gz     |
| ERR3333390 | R        | ERR3333390_R1  | transcriptome| ERR3333390_q30l50_R1.fastq.gz     |

Is just an inventory of dataset filenames used for each Group defined in the file of files provided to Shake by the user
