# Shake - Shannon diversity on K-mers and gene Expression

## Requirements

**smudgeplot v0.2.5** : We recommend to use conda (`conda install smudgeplot==0.2.5`) to install this module.  

kmc https://github.com/refresh-bio/KMC. You will need to add kmc binaries to PATH variable for example using `export PATH=$PATH:kmc/install/bin/`

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

The third column serves to label what kind of data is; you can put whaterver you like but if you use *** transcriptome *** hetkmers will be calculated and kmers and their values will be grouped using this. If you put *** expression ***, files will be treated as expression tables.

Fourth and fifth column are the minimun and the maximum values cutoff in order to consider a kmer or not in the analysis.

Sixth column is the filepaths column, you can separate different files with a comma.

This pipeline can be run with `omics_diversity_pipeline.py`and it accepts the following arguments:

  `--input_file INPUT_FILE, -i INPUT_FILE` (Required): the path to the table previously described
  
  `--output_dir OUTPUT_DIR, -o OUTPUT_DIR`(Required): the path were interemediate files and results will be stored. It will also check in this ouput for steps already done and will skip these steps.
  
  `--ram_usage RAM_USAGE, -r RAM_USAGE` (Optional, 1 by default): The max memory usage of kmc.
  
  `--num_threads NUM_THREADS, -t NUM_THREADS`(Optional 1 by default): Number of threads used by kmc.
  
  `--kmer_size KMER_SIZE, -k KMER_SIZE` (Optional, 21 by default): the kmer size used by kmc.
  
  `--merge_universe, -m` (Optional, False by default): This switchs whether if the kmer universe size should be calculated from all items from a group (True) or only for each sub group.
  
  `--presence, -p` (Optional, False by default) switchs change diversity calculation to presence/absence, counting all kmer count values to 1 if is present. For expression tables, a gene will throw a value of 1 if a gene has a TPM of 1 or more or 0 if its less than 1. 

  ## Output

  ### Diversity results

  Shake will generate a table called `results.tsv` inside the directory specified with the `--output_dir` argument. This table has the following format:

| Group      | Subgroup | Rep              | Kind         | Subgroup_Universe_Size | Diversity_log2 | Specifity_log2 | Diversity_log10 | Specifity_log10 | Kolmogorov         |
|------------|----------|------------------|--------------|------------------------|----------------|----------------|------------------|------------------|---------------------|
| ERR3333388 | R        | ERR3333388_R1    | transcriptome| 40072639               | 22.7620354967  | 2.4940786897   | 6.8520554469     | 0.7507924971     | 0.0636632271        |
| ERR3333389 | R        | ERR3333389_R1    | transcriptome| 38691164               | 22.9552288052  | 2.2502719911   | 6.9102124277     | 0.6773993677     | 0.0648244898        |
| ERR3333390 | R        | ERR3333390_R1    | transcriptome| 39507337               | 22.9614172079  | 2.2742000608   | 6.9120753225     | 0.6846024344     | 0.0653194914        |
| ERR3333391 | R        | ERR3333391_R1    | transcriptome| 40188471               | 22.5576212642  | 2.7026570900   | 6.7905206314     | 0.8135808521     | 0.0647497717        |
| ERR3333392 | R        | ERR3333392_R1    | transcriptome| 36116326               | 21.8123101535  | 3.2938376502   | 6.5661596309     | 0.9915439335     | 0.0630993513        |


  `Group , Subgroup, Rep and Kind`: have the values provided by the user in the file of files.
  
  `Subgroup_Universe_Size`: is the number of unique K-mers (universe size) found in the dataset.
  
  `Diversity_log2 and Diversity_log10`: is the Shannon Diversity Index score calculated to for each dataset, using log2 and log10.
  
  `Specificty_log2 and Specifity_log10`:is the K-mer specifity of each dataset. A specifity value of 0 means that each kind of K-mer is equally represented in a given dataset and higher values means that some K-mers are overrepresented. Specificty is calculated for log2 and log10. It would be 0 or near 0 `--presence` option is used
  
  `Kolmogorov`: is Kolmogorov complexity measurement. Is an alternative measurement for Shannon Diversity Index that reduces the impact of K-mers exclusive to unique datasets from the same group.

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
