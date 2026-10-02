#!/usr/bin/env python
"""Generate the file of files (fof) that SHaKE.py takes with --input_file.

Every line of the fof is tab separated:

    group  subgroup  kind  lowerbound  upperbound  file1[,file2]

This script scans a directory of reads and writes one line per sample.
Paired reads (two files per sample, e.g. sample_1.fastq.gz / sample_2.fastq.gz,
sample_R1.fq.gz / sample_R2.fq.gz, sample_R1_001.fastq.gz / sample_R2_001.fastq.gz)
are put together on the same line, comma separated, read 1 first. Files
without a mate are written alone (single-end reads, fasta assemblies, bam).

Usage:
    python make_fof.py -i /path/to/reads -o fof.tsv --subgroup R --kind transcriptome
    python make_fof.py -i /path/to/reads --group onegroup --subgroup R --kind genome
"""
import argparse
import re
import sys

from pathlib import Path

SEQ_EXTENSIONS = r"(?:fastq|fq|fasta|fa|fna|bam)(?:\.gz)?"
# <sample><sep><1|2><optional _001><extension>; sep is one of _ . - optionally followed by R
PAIRED_RE = re.compile(r"^(?P<sample>.+?)(?P<sep>[_.-]R?)(?P<read>[12])(?P<tail>_001)?"
                       r"\.(?P<ext>" + SEQ_EXTENSIONS + r")$", re.IGNORECASE)
SEQ_RE = re.compile(r"^(?P<sample>.+)\.(?P<ext>" + SEQ_EXTENSIONS + r")$", re.IGNORECASE)


def parse_arguments():
    desc = "Generate the file of files (fof) used as input by SHaKE.py"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--input_dir", "-i", type=str, required=True,
                        help="(Required) directory with the reads")
    parser.add_argument("--output", "-o", type=str, default="",
                        help="(Optional) fof path to write. Printed to stdout by default")
    parser.add_argument("--group", "-g", type=str, default="",
                        help="(Optional) value for the group column, shared by all "
                             "samples. By default each sample is its own group")
    parser.add_argument("--subgroup", "-s", type=str, default="R",
                        help="(Optional) value for the subgroup column. R by default")
    parser.add_argument("--kind", "-k", type=str, default="transcriptome",
                        help="(Optional) data kind column (transcriptome, genome, "
                             "expression...). transcriptome by default")
    parser.add_argument("--lowerbound", "-l", type=int, default=0,
                        help="(Optional) minimum kmer count cutoff. 0 by default")
    parser.add_argument("--upperbound", "-u", type=int, default=9999999999,
                        help="(Optional) maximum kmer count cutoff. 9999999999 by default")
    parser.add_argument("--recursive", "-r", action="store_true", default=False,
                        help="(Optional) also search subdirectories. False by default")
    parser.add_argument("--relative", action="store_true", default=False,
                        help="(Optional) write paths as given instead of absolute ones")
    parser.add_argument("--glob", type=str, default="",
                        help="(Optional) only use files whose name matches this "
                             "pattern, e.g. '*_q30l50_*'")
    return parser.parse_args()


def find_reads(input_dir, recursive=False, pattern=""):
    paths = input_dir.rglob("*") if recursive else input_dir.glob("*")
    for path in sorted(paths):
        if not path.is_file() or not SEQ_RE.match(path.name):
            continue
        if pattern and not path.match(pattern):
            continue
        yield path


def group_samples(paths):
    """Returns a list of (sample_name, [files]) with the mates of paired reads together"""
    pairs = {}
    singles = []
    for path in paths:
        match = PAIRED_RE.match(path.name)
        if not match:
            singles.append((SEQ_RE.match(path.name).group("sample"), path))
            continue
        # files are mates only if everything but the 1/2 is identical
        key = (path.parent, match.group("sample"), match.group("sep").upper(),
               match.group("tail") or "", match.group("ext").lower())
        pairs.setdefault(key, {})[match.group("read")] = path

    samples = []
    for key, reads in pairs.items():
        sample = key[1]
        if len(reads) == 2:
            samples.append((sample, [reads["1"], reads["2"]]))
        else:
            read, path = next(iter(reads.items()))
            print("WARNING: {} looks like read {} but its mate was not found, "
                  "using it as single file".format(path, read), file=sys.stderr)
            samples.append((sample, [path]))
    samples += [(sample, [path]) for sample, path in singles]
    return sorted(samples, key=lambda sample: (sample[0], str(sample[1][0])))


def main():
    args = parse_arguments()
    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        sys.exit("ERROR: {} is not a directory".format(input_dir))

    samples = group_samples(find_reads(input_dir, args.recursive, args.glob))
    if not samples:
        sys.exit("ERROR: no fasta/fastq/bam files found in {}".format(input_dir))

    lines = []
    for sample, files in samples:
        group = args.group if args.group else sample
        paths = [str(path if args.relative else path.resolve()) for path in files]
        lines.append("\t".join([group, args.subgroup, args.kind,
                                str(args.lowerbound), str(args.upperbound),
                                ",".join(paths)]))

    text = "\n".join(lines) + "\n"
    if args.output:
        Path(args.output).write_text(text)
        n_paired = sum(1 for _, files in samples if len(files) == 2)
        print("Wrote {} samples ({} paired, {} single) to {}".format(
            len(samples), n_paired, len(samples) - n_paired, args.output), file=sys.stderr)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
