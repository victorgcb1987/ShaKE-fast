#!/usr/bin/env python
"""Generate a synthetic transcriptome dataset with known het-kmer scenarios,
to validate SHaKE's native het-kmer detection/merging (src/kmerdata.py) end
to end through the real kmc/kmc_tools pipeline - not just the internal unit
checks.

Produces four independent loci in one FASTA file, each built so its expected
outcome is known in advance:

  balanced_snp    one SNP, ~equal allele depth (20 vs 20)
                  -> should MERGE in both --pooled and default mode
                     (component size 2, coverage ratio ~1.0)

  imbalanced_snp  one SNP, lopsided allele depth (20 vs 1)
                  -> should merge ONLY with --pooled
                     (component size 2 either way, but ratio 0.05 fails the
                      non-pooled coverage-ratio check)

  repeat_hub      one "hub" k-mer with 7 single-substitution "leaf" neighbors,
                  mimicking a repeat/paralog family
                  -> should NEVER merge in either mode
                     (component size 8 exceeds both the pooled cap of 6 and
                      the non-pooled cap of 2)

  invariant       one plain locus with no variants at all
                  -> negative control: no het-kmer edges at all, must be left
                     completely untouched

Also writes a matching file-of-files ready to hand to SHaKE.py, and prints
the exact numbers you should see in the run's log and results.tsv.

Usage:
    python validation/make_synthetic_hetkmer_dataset.py -o /tmp/hetkmer_test
    python SHaKE.py -i /tmp/hetkmer_test/file_of_files.tsv -o /tmp/hetkmer_test/run_default
    python SHaKE.py -i /tmp/hetkmer_test/file_of_files.tsv -o /tmp/hetkmer_test/run_pooled --pooled
"""
import argparse
import random

from pathlib import Path

K = 21
N_LEAVES = 7


def random_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def substitute(seq, pos, new_base):
    return seq[:pos] + new_base + seq[pos + 1:]


def other_base(base, rng):
    return rng.choice([b for b in "ACGT" if b != base])


def write_reads(fhand, seq, count, name_prefix, start_index):
    for i in range(count):
        fhand.write(">{}_{}\n{}\n".format(name_prefix, start_index + i, seq))
    return start_index + count


def build_fasta(fasta_path, rng):
    with open(fasta_path, "w") as fhand:
        idx = 0

        ref = random_seq(K, rng)
        alt = substitute(ref, K // 2, other_base(ref[K // 2], rng))
        idx = write_reads(fhand, ref, 20, "balanced_ref", idx)
        idx = write_reads(fhand, alt, 20, "balanced_alt", idx)

        ref2 = random_seq(K, rng)
        alt2 = substitute(ref2, K // 2, other_base(ref2[K // 2], rng))
        idx = write_reads(fhand, ref2, 20, "imbalanced_ref", idx)
        idx = write_reads(fhand, alt2, 1, "imbalanced_alt", idx)

        hub = random_seq(K, rng)
        idx = write_reads(fhand, hub, 10, "hub", idx)
        leaf_positions = [pos for pos in range(0, K, 3)][:N_LEAVES]
        for leaf_i, pos in enumerate(leaf_positions):
            leaf = substitute(hub, pos, other_base(hub[pos], rng))
            idx = write_reads(fhand, leaf, 10, "leaf{}".format(leaf_i), idx)

        invariant = random_seq(K, rng)
        idx = write_reads(fhand, invariant, 15, "invariant", idx)
    return idx


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output_dir", "-o", required=True)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / "synthetic_hetkmer_test.fasta"
    fof_path = out_dir / "file_of_files.tsv"

    rng = random.Random(args.seed)
    n_reads = build_fasta(fasta_path, rng)

    with open(fof_path, "w") as fhand:
        fhand.write("hetkmer_test\tsample1\ttranscriptome\t1\t9999999999\t{}\n".format(fasta_path))

    n_leaf_members = 1 + N_LEAVES  # hub + leaves
    raw_universe = 2 + 2 + n_leaf_members + 1  # balanced + imbalanced + repeat_hub + invariant

    print("Wrote {} reads to {}".format(n_reads, fasta_path))
    print("Wrote file-of-files to {}".format(fof_path))
    print()
    print("Run SHaKE.py twice against this file-of-files (default, then --pooled) and compare:")
    print()
    print("Expected STEP 4 log line (grep '#SUCCESS.*components_found'):")
    print("  components_found=3  (balanced_snp, imbalanced_snp, repeat_hub each form one")
    print("                       component; invariant has no edges so contributes none)")
    print("  default   : components_merged=1 components_rejected=2")
    print("              (balanced_snp merges; imbalanced_snp fails the coverage-ratio")
    print("               check; repeat_hub's size {} exceeds the non-pooled cap of 2)".format(n_leaf_members))
    print("  --pooled  : components_merged=2 components_rejected=1")
    print("              (balanced_snp AND imbalanced_snp merge - no ratio check applies;")
    print("               repeat_hub's size {} still exceeds the pooled cap of 6)".format(n_leaf_members))
    print()
    print("Expected Subgroup_Universe_Size in results.tsv:")
    print("  raw (no merging)       : {}".format(raw_universe))
    print("  default (1 merge)      : {}".format(raw_universe - 1))
    print("  --pooled (2 merges)    : {}".format(raw_universe - 2))


if __name__ == "__main__":
    main()
