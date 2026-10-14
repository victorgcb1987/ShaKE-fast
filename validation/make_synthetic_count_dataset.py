#!/usr/bin/env python
"""Generate a small synthetic transcriptome dataset whose k-mer count histograms are
known exactly, to validate the count-based complexity columns of results.tsv
(Spectrum_entropy_universe, Occupancy_entropy, Spectrum_entropy, Spectrum_entropy_norm).

Three samples are built from random "transcripts" (random 120-mers, whose 21-mers
are unique and do not overlap between transcripts), each transcript sequenced as
N identical reads, so every k-mer of a transcript has count N:

  even      transcripts T00-T39, 10 reads each            (even counts)
  skewed    T00-T19 (shared with `even`) plus T40-T59 new;
            T40 has 200 reads, every other transcript 5   (one dominant transcript)
  small     T00-T09 only, 10 reads each                   (subset of `even`)

Two files of files are written:
  fof_shared.tsv    the three samples in ONE group/subgroup, so the universe is the
                    union of their k-mers (the zeros are padded to that size)
  fof_separate.tsv  each sample in its own group (universe = the sample's own k-mers,
                    no zeros: what happens when every sample is its own group)

The script prints the values each column must have, computed independently from the
count histograms (direct sum, not through src/kmerdata.py):

    Spectrum_entropy_universe     = sum_c N_c * log2(U / N_c) / U      (c = 0 included, N_0 = U - n)
    Occupancy_entropy  = H2(n / U)
    Spectrum_entropy        = H(count values among the present k-mers)
    Spectrum_entropy_norm= Spectrum_entropy / Hmax(mean count - minimum count)

With --check the expected values are compared with a results.tsv (max absolute
difference per column, tolerance 1e-9); the row "group_R_merged" (the pooled database of the
subgroup, only present in the shared layout) is not checked:
    python validation/make_synthetic_count_dataset.py -o /tmp/count_test \
        --check /tmp/count_test/run_shared/results.tsv --layout shared

Usage:
    python validation/make_synthetic_count_dataset.py -o /tmp/count_test
    python SHaKE.py -i /tmp/count_test/fof_shared.tsv -o /tmp/count_test/run_shared
    python SHaKE.py -i /tmp/count_test/fof_separate.tsv -o /tmp/count_test/run_separate
    python SHaKE.py -i /tmp/count_test/fof_shared.tsv -o /tmp/count_test/run_shared_seq --sequential
"""
import argparse
import math
import random

from collections import Counter
from pathlib import Path

K = 21
LENGTH = 120


def random_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def kmers(seq):
    return [seq[i:i + K] for i in range(len(seq) - K + 1)]


def h2(p):
    return 0.0 if p <= 0 or p >= 1 else -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def hmax(mu):
    # maximum entropy of a distribution on {0, 1, 2, ...} with mean mu (geometric)
    return 0.0 if mu <= 0 else (mu + 1) * math.log2(mu + 1) - mu * math.log2(mu)


def expected(counts, universe):
    """counts: dict kmer -> count (all > 0). universe: number of elements."""
    hist = Counter(counts.values())
    n = sum(hist.values())
    zeros = universe - n
    direct = sum(f * math.log2(universe / f) for f in hist.values())
    if zeros > 0:
        direct += zeros * math.log2(universe / zeros)
    h_c = abs(-sum(f / n * math.log2(f / n) for f in hist.values()))
    mean = sum(c * f for c, f in hist.items()) / n
    norm = h_c / hmax(mean - min(hist)) if hmax(mean - min(hist)) > 0 else float("nan")
    return {"n": n, "U": universe,
            "Spectrum_entropy_universe": direct / universe,
            "Occupancy_entropy": h2(n / universe),
            "Spectrum_entropy": h_c,
            "Spectrum_entropy_norm": norm}


COLUMNS = {"Spectrum_entropy_universe": "Spectrum_entropy_universe", "Occupancy_entropy": "Occupancy_entropy",
           "Spectrum_entropy": "Spectrum_entropy", "Spectrum_entropy_norm": "Spectrum_entropy_norm"}


def check(results_path, layout, counts, union_size):
    import csv
    rows = {row["Rep"]: row for row in csv.DictReader(open(results_path), delimiter="\t")}
    worst = 0.0
    for i, (name, c) in enumerate(counts.items(), start=1):
        rep = "group_R{}".format(i) if layout == "shared" else "{}_R1".format(name)
        e = expected(c, union_size if layout == "shared" else len(c))
        for column in COLUMNS:
            got, want = float(rows[rep][column]), e[column]
            if math.isnan(want):
                ok = math.isnan(got)
                diff = 0.0 if ok else float("inf")
            else:
                diff = abs(got - want)
                ok = diff < 1e-9
            worst = max(worst, diff)
            print("  {:<8} {:<22} expected {:<12.9f} got {:<12.9f} {}".format(name, column, want, got, "ok" if ok else "MISMATCH"))
    print("max absolute difference: {:.2e}".format(worst))
    raise SystemExit(0 if worst < 1e-9 else 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--check", help="results.tsv of a SHaKE.py run to compare with the expected values")
    parser.add_argument("--layout", choices=["shared", "separate"], default="shared",
                        help="which file of files the --check results.tsv comes from")
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    transcripts = ["T{:02d}".format(i) for i in range(60)]
    seqs = {t: random_seq(LENGTH, rng) for t in transcripts}
    samples = {
        "even": {"T{:02d}".format(i): 10 for i in range(40)},
        "skewed": dict([("T{:02d}".format(i), 5) for i in range(20)] +
                       [("T{:02d}".format(i), 200 if i == 40 else 5) for i in range(40, 60)]),
        "small": {"T{:02d}".format(i): 10 for i in range(10)},
    }
    counts = {}
    for name, reads in samples.items():
        with open(out / "{}.fasta".format(name), "w") as fhand:
            for t, n_reads in reads.items():
                for j in range(n_reads):
                    fhand.write(">{}_{}\n{}\n".format(t, j, seqs[t]))
        counts[name] = {}
        for t, n_reads in reads.items():
            for km in kmers(seqs[t]):
                counts[name][km] = counts[name].get(km, 0) + n_reads
    with open(out / "fof_shared.tsv", "w") as shared, open(out / "fof_separate.tsv", "w") as separate:
        for name in samples:
            path = (out / "{}.fasta".format(name)).resolve()
            shared.write("group\tR\ttranscriptome\t1\t9999999999\t{}\n".format(path))
            separate.write("{}\tR\ttranscriptome\t1\t9999999999\t{}\n".format(name, path))
    union = set()
    for c in counts.values():
        union |= set(c)
    if args.check:
        check(args.check, args.layout, counts, len(union))
        return
    print("expected values (k = {}; union of k-mers = {})".format(K, len(union)))
    for layout, universe in (("shared (fof_shared.tsv)", lambda c: len(union)),
                             ("separate (fof_separate.tsv)", lambda c: len(c))):
        print("\n{}".format(layout))
        for name, c in counts.items():
            e = expected(c, universe(c))
            print("  {:<7} n={:<5} U={:<5} Spectrum_entropy_universe={:.6f} Occupancy_entropy={:.6f} "
                  "Spectrum_entropy={:.6f} Spectrum_entropy_norm={:.6f}".format(
                      name, e["n"], e["U"], e["Spectrum_entropy_universe"], e["Occupancy_entropy"],
                      e["Spectrum_entropy"], e["Spectrum_entropy_norm"]))


if __name__ == "__main__":
    main()
