# Methodology notes

This document explains **why** ShaKE's pipeline looks the way it does, from both
a biological and a software-engineering angle. `CHANGELOG.md` records *what*
changed line-by-line; `PERFORMANCE.md` covers an earlier, narrower round of
speed fixes. This document is for the bigger question: why does the code make
these particular scientific and design choices, so a future reader (including
future-you) doesn't have to reverse-engineer the reasoning from a diff.

## 1. What ShaKE actually computes

ShaKE treats a sample's k-mer counts (or, for expression data, a gene's TPM
values) the way ecology treats species counts in a quadrat: each distinct
k-mer/gene is a "species," its count is its abundance, and the pipeline
reports:

- **Shannon diversity** (`Diversity_log2`/`Diversity_log10`) - how evenly
  spread the abundance is across all the distinct k-mers/genes observed.
  Low diversity means a few k-mers dominate; high diversity means abundance is
  spread evenly across many.
- **Specificity** (`Specifity_log2`/`Specifity_log10`) - `log(universe_size) -
  diversity`. This is "how far from maximally even" the sample is, scaled by
  how big the universe of possible k-mers/genes is.
- **Kolmogorov complexity ratio** - an alternative, non-entropy-based
  diversity proxy: encode each count as a fixed-width bit string, compress the
  result, and report `compressed_size / uncompressed_size`. A sample whose
  counts are highly repetitive/predictable compresses well (low ratio); one
  whose counts vary a lot compresses poorly (high ratio). It's included
  because it reacts differently than Shannon entropy to k-mers that are
  exclusive to one sample in a group, which is useful when comparing samples
  within the same group.

Every one of these is now computed twice per sample - once from raw counts,
once from presence/absence (every observed k-mer/gene counted as 1 regardless
of its actual abundance) - and reported side by side as `..._presence`
columns. See §3 for why.

## 2. Structural refactor (informatics)

The pipeline used to live entirely inside `SHaKE.py`'s `main()`, as one
~270-line sequential function sharing two large mutable dicts across five
"STEP" blocks. This was reorganized into named stage functions in
`src/pipeline.py` (`build_databases`, `build_histograms`, `dump_counts`,
`merge_hetkmers`, `compute_universe_sizes`, `compute_estimators`,
`compute_expression_estimators`, `write_outputs`), with `SHaKE.py` reduced to
argument parsing plus a call into `run_pipeline()`. No behavior was
intentionally changed by this reorganization - it's purely so each step can be
read, tested, and modified independently. See `CHANGELOG.md` for the exact
before/after.

Alongside the reorganization, four small pre-existing bugs were found and
fixed (each was surfaced and confirmed individually rather than silently
folded in - see `CHANGELOG.md` for details):
- `expression`-kind rows in `results.tsv` were reporting a stale, unrelated
  sample's numbers instead of their own.
- `file_manifiest.tsv`'s header was malformed (missing a tab, stray quote).
- `main()` called the input-parsing function twice, wastefully re-running
  BAM→FASTA conversion and leaking a file handle.
- Several "already done, skip on rerun" checks reported inconsistently (some
  logged nothing at all).

One bug was found and *deliberately left alone*: the k-mer dump step reads its
count-cutoff bounds from the wrong dictionary level, so it silently falls back
to defaults instead of the input file's actual bounds. It's not fixed because
the k-mer *counting* step upstream already enforces the real bounds correctly,
so this is likely harmless double-filtering rather than a live bug - but
changing it would alter which k-mers appear in the dump, so it's flagged
rather than touched without being asked.

## 3. Reporting both count-based and presence/absence diversity

**The `--presence` flag was removed; both variants are now always computed.**
Rationale: count-based and presence/absence diversity answer genuinely
different biological questions - "how evenly is abundance spread" versus "how
many distinct things are present at all, ignoring how much of each." Making
the user pick one upfront via a flag means re-running the whole pipeline to
get the other view. Since both are cheap to compute once you already have the
counts in hand, there's no real cost to always reporting both, and it removes
a decision the user shouldn't have to make in advance.

## 4. Native het-kmer detection for pooled resequencing

This is the part of the pipeline most directly shaped by the specific biology
of the data it's run on, so it gets the longest explanation.

### 4.1 The biological problem

A resequencing project that pools multiple individuals into one sample means
a single genomic/transcript locus can contribute **more than one distinct
k-mer string** to the count table, purely because a SNP segregating across
the pooled individuals falls inside that k-mer's window. Without correction,
a k-mer-based diversity index treats each of those SNP-variant k-mers as a
separate "species," inflating apparent diversity/richness for reasons that
have nothing to do with the actual signal of interest (e.g. tissue
specificity, expression breadth). The fix is to detect k-mers that differ by
exactly one substitution and are otherwise the same, and collapse them back
into one unit (summing their counts) before computing diversity - this is
what "het-kmer" detection and merging does.

### 4.2 Why the previous approach (`smudgeplot.py hetkmers`) didn't fit

Smudgeplot's het-kmer detection was built for a different question: calling
heterozygous sites in **one diploid individual**. Its model expects a
candidate het-kmer pair to have roughly balanced coverage (~50/50), because
that's what two alleles of one het site look like in a single genome
sequenced at uniform depth. In a pool of many individuals:
- allele frequency at a segregating site can be anything, not just ~50/50 -
  a rare variant present in one of many pooled individuals might show a
  95/5 split, and smudgeplot's balance assumption would reject or mis-score
  it;
- more than two alleles can legitimately co-occur at one site across a pool,
  which a strictly-pairwise heterozygous model doesn't represent well.

So smudgeplot wasn't just an external dependency to remove for its own sake -
it was actively the wrong statistical model for this data.

### 4.3 The replacement: native Hamming-distance-1 detection

Implemented in `src/kmerdata.py`. The core idea: two k-mers that differ by
exactly one substitution are Hamming-distance-1 neighbors. For each observed
k-mer, generate its `3k` possible single-substitution variants (3 alternate
bases at each of `k` positions) and check whether any of them is *also*
present in the sample's k-mer set. Every hit becomes an edge in a graph;
connected components of that graph are groups of k-mers hypothesized to
originate from the same locus. This needs no assumption about coverage
balance at all - it's purely "are these two k-mer strings one substitution
apart," which is exactly the definition of a candidate SNP pair.

Practically, this is done with k-mers 2-bit-encoded as integers (A/C/G/T →
0/1/2/3, packed MSB-first into a 64-bit int, which fits any k ≤ 32) so that
neighbor generation is XOR against a precomputed bitmask and membership
testing is a vectorized numpy binary search, instead of doing 3k string
substitutions and dict lookups per k-mer in a Python loop. At the scale a
pooled dataset can reach (tens to hundreds of millions of k-mers), this
difference is the difference between "runs" and "doesn't."

One correctness detail worth remembering: KMC counts a k-mer and its reverse
complement as the same entity, keeping whichever is lexicographically smaller
("canonical form"). Since the 2-bit encoding preserves string ordering, this
is reproduced exactly by `min(encoded, revcomp(encoded))` on the integers -
every substitution candidate is re-canonicalized this way before the
membership check, or real SNP pairs would be missed whenever a substitution
happens to flip which strand is canonical.

### 4.4 Why k=21 mostly - but not entirely - solves the "is this really the
same locus" question

A natural objection: without a coverage-balance signal, how do you know a
Hamming-1 pair reflects the *same genomic origin*, rather than two unrelated
k-mers that just happen to be one substitution apart? The answer has two
parts:

**The k-mer space is enormous relative to any real dataset.** At k=21, there
are 4²¹ ≈ 4.4×10¹² possible k-mer sequences, and each k-mer has only `3×21 =
63` possible single-substitution neighbors. For a specific k-mer, the odds
that some *unrelated* k-mer elsewhere in a transcriptome of size N happens to
be one of those 63 sequences is roughly `63×N / 4^21`. For N≈5×10⁷ (a
sizeable pooled transcriptome) that's about 0.07% per k-mer - this is
precisely why k≈21 is the conventional choice for this class of analysis
(GenomeScope/Merqury/smudgeplot all lean on the same argument): long enough
to make random collisions between unrelated loci rare, short enough to still
tolerate one SNP without losing all overlap.

**But "rare per k-mer" isn't "rare overall" at pool scale**, and there's a
second, non-random confound k=21 doesn't touch at all. Working the same
number through for a whole dataset (`N × 63 × N / 4^21 / 2` for N≈5×10⁷)
gives roughly 18,000 expected coincidental pairs - small as a fraction of
the total, but not zero. More importantly, **repeats and paralogs** (gene
families, recent duplications, transposable elements) produce genuinely
distinct genomic loci that are nearly identical in sequence for reasons that
have nothing to do with allelic variation - and coverage balance was never
protecting against this either, since a paralog can easily have coverage
similar to its sibling copy.

### 4.5 The mitigation: component-size capping, switched by `--pooled`

Since repeats/paralogs are the dominant real risk (not random collision), the
chosen mitigation is topological rather than coverage-based: **a true SNP
site should only ever produce a handful of alleles**, so a connected
component in the Hamming-1 graph larger than that ceiling is treated as a
repeat-family signature and left un-merged, regardless of mode.

- **`--pooled`** (the primary intended use case - multiple individuals pooled
  into one resequencing sample): component size capped at **6** (generous
  enough for several co-segregating alleles across a pool), **no**
  coverage-ratio requirement (allele frequency in a pool can legitimately be
  anything).
- **default / not `--pooled`** (a single individual): component size capped
  at exactly **2** (a single diploid genome is strictly biallelic at any
  site, barring CNVs), **plus** a coverage-ratio check
  (`min(count)/max(count) ≥ 0.2`) - this is the regime where coverage balance
  is actually a valid signal again (it's the model smudgeplot itself
  assumed), so it's kept as an extra check rather than discarded. A pair
  failing the ratio check is presumed more likely to be a sequencing artifact
  than a real het site and is left un-merged.

Both the cap values and the ratio threshold are plain constants in
`src/kmerdata.py` (`POOLED_MAX_COMPONENT_SIZE`, `NONPOOLED_MAX_COMPONENT_SIZE`,
`NONPOOLED_MIN_COVERAGE_RATIO`), not exposed as further CLI flags - the
intent was to keep the user-facing surface to the one `--pooled` switch,
since these are judgment calls about acceptable false-positive/false-negative
trade-offs rather than something to tune per run.

### 4.6 What this still doesn't solve

Two things worth remembering if diversity numbers ever look off on real data:

- **Sequencing errors near the count cutoff.** Every single-base sequencing
  error manufactures a Hamming-1 neighbor of the true k-mer, but at low
  count - the pipeline's existing `lowerbound`/`upperbound` cutoffs (applied
  at k-mer counting time) filter most of this out before het-kmer detection
  ever runs. The one case that stays genuinely hard is a *rare real* pooled
  allele sitting close to that same low-count cutoff, which can't be cleanly
  distinguished from a sequencing error by count alone.
- **Multiple SNPs inside one k-mer window**, or a locus with more genuine
  variation than the component-size cap allows, won't be recognized as one
  site by this method. The more rigorous alternative, if this ever becomes a
  real limitation, is de Bruijn "bubble" detection (the approach tools like
  DiscoSNP++/KisSNP use): instead of trusting an isolated Hamming-1 pair,
  require that the two candidate k-mers sit on two paths through the de
  Bruijn graph that diverge at one node and reconverge k steps later - i.e.
  corroborate the SNP with agreement on the *flanking* sequence, not just the
  substituted base. That needs k-mer adjacency information (which base
  extends a k-mer to k+1), not just a count table, so it's a meaningfully
  bigger piece of work than what's implemented here - noted as a possible
  future direction, not attempted.

## 5. Disk footprint and streaming redesign (informatics)

### 5.1 The problem

Every sample used to go through: `kmc_tools dump` writing a full text file of
every k-mer and its count → read back separately by three different pipeline
stages (universe size, het-kmer merge, Shannon estimator) → a further
`.binary`/`.presence.binary` file (one line per k-mer, encoded as a
fixed-width bit string) plus its gzip'd copy, written *twice* per sample once
the presence/absence feature (§3) doubled that step. None of these
intermediate files are useful outputs in their own right - they exist only to
get one number (a count, an entropy value, a compression ratio) out the other
side.

### 5.2 The fix

`kmc_tools ... dump -s /dev/stdout` streams the dump output through a
subprocess pipe straight into Python (confirmed empirically to leave no file
behind, both via a shell pipe and via `subprocess.Popen(..., stdout=PIPE)`),
building one compact in-memory table (`KmerTable` in `src/kmerdata.py`:
sorted numpy arrays of 2-bit-encoded k-mer integers + counts). That single
table is then reused for universe size (its length), Shannon entropy (a
vectorized numpy computation over the counts array), het-kmer detection (§4),
and the Kolmogorov ratio - the fixed-width bit-string encoding is built as an
in-memory byte buffer and compressed with `zlib.compress()` directly, so
nothing is ever written to disk for that step either, for genomic,
transcriptome, or expression samples.

### 5.3 The trade-off, taken on deliberately

The existing "if the output already exists, skip this step" resumability
depended on those intermediate files being there to check. Once they're gone,
there's nothing left to check, so the dump/het-kmer-merge/estimator/Kolmogorov
steps now **always recompute** on every run. This was accepted explicitly
rather than solved with a different checkpointing mechanism, because:
- the expensive step - kmc's own database build (`.kmc_pre`/`.kmc_suf`,
  which is where the real, size-of-the-raw-sequencing-data cost lives) still
  gets skipped on a rerun, unaffected by any of this;
- everything downstream of that (dump → estimators) is now fast enough,
  running from an already-built kmc database, that recomputing it every time
  is a reasonable price for not accumulating disk usage or maintaining a
  second checkpointing scheme for cheap steps.

## 6. Known, deliberately-not-fixed loose ends

Kept here so they don't get "rediscovered" as surprises later:

- **Dump-step cutoff bug** (§2) - reads bounds from the wrong dict level,
  falls back to defaults. Left alone because kmc counting upstream already
  enforces the real bounds.
- **Kolmogorov compression-level inconsistency** - the genomic/transcriptome
  path uses zlib level 6 (matching the old `gzip -c` shell command's
  default); the expression path uses level 9 (matching the old Python
  `gzip.open()` default). This predates the streaming redesign - it was
  simply two different pieces of code that happened to pick different
  defaults - and was preserved rather than unified, so that switching to the
  in-memory implementation didn't also silently change expression samples'
  Kolmogorov numbers relative to prior runs. Worth unifying deliberately at
  some point if genomic and expression Kolmogorov values are ever compared
  directly against each other.
- **Old file-based Shannon/Kolmogorov functions** (`src/kmer.py:
  calculate_sample_shannon_estimators`, `src/kolmogorov.py:
  calculate_kolmogorov_estimator` and its helpers) are no longer called by
  the pipeline but were kept in place specifically to make old-vs-new parity
  testing possible. They're dead code from the pipeline's point of view -
  reasonable to remove in a future cleanup pass once nobody needs to diff
  against them anymore.

## 7. Quick map of where things live

| File | Responsibility |
|---|---|
| `SHaKE.py` | CLI argument parsing only; calls `src.pipeline.run_pipeline` |
| `src/pipeline.py` | Stage orchestration - one function per pipeline step |
| `src/kmc.py` | Wraps external `kmc`/`kmc_tools` calls (database build, histogram, streaming dump command) |
| `src/kmerdata.py` | In-memory k-mer table, encoding, native het-kmer detection, Shannon entropy, Kolmogorov ratio - the numerical core |
| `src/kolmogorov.py` | Expression-path Kolmogorov ratio; older file-based functions kept for parity testing only |
| `src/kmer.py` | Older file-based Shannon estimator, kept for parity testing only |
| `src/expression.py` | Expression-table (TPM) Shannon estimator |
| `src/utils.py` | Shared helpers: union-find, BAM/FASTA handling, status logging |
| `legacy_scripts/` | Predecessor standalone scripts; untouched, still using the older `src/kmc.py`/`src/utils.py` functions directly |
