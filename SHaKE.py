#!/usr/bin/env python


import argparse
import sys

from datetime import datetime
from pathlib import Path

from src.utils import check_run, sequence_kind, convert_bam_to_fasta, filter_bam_file
from src.pipeline import run_pipeline


def parse_arguments():
    desc = "Pipeline to run all steps required"
    parser = argparse.ArgumentParser(description=desc)


    help_input = '''(Required) File of files including kind
                    Name1  R1   genomic  KmerLength lower_bound upper_bound file1,file2
                    '''
    parser.add_argument("--input_file", "-i", type=str,
                        help=help_input, required=True)

    help_output = "(Required) output dir"
    parser.add_argument("--output_dir", "-o", type=str,
                        help=help_output, required=True)

    help_ram = "(Optional) Max RAM usage. 6GB by default"
    parser.add_argument("--ram_usage", "-r", type=int,
                        help=help_ram, default=6)
    help_threads = "(Optional) Number of threads. 1 by default"
    parser.add_argument("--num_threads", "-t", type=int,
                        help=help_threads, default=1)
    help_kmer_size = "(Optional) Kmer size. 21 by default"
    parser.add_argument("--kmer_size", "-k", type=int,
                        help=help_kmer_size, default=21)
    help_universe = "(Optional) merge universe within a group ignoring subrgroups. False by default"
    parser.add_argument("--merge_universe", "-m",
                        help=help_universe, default=False,
                        action="store_true")
    help_exclude = "(optional) Bed file with regions to exclude"
    parser.add_argument("--exclude", "-e",
                        help=help_exclude, type=str,
                        default="")


    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    inputs = {}
    input_sequence = Path(parser.input_file)
    logdate = "Kmer_counting_"+datetime.now().strftime("%d_%m_%Y-%H_%M_%S") + ".log"
    if not Path(parser.output_dir).exists():
        Path(parser.output_dir).mkdir(exist_ok=True, parents=True)
    log_fname = Path(parser.output_dir) / logdate
    log_fhand = open(log_fname, "w")
    msg = "#Command Used: "+ " ".join(sys.argv)+"\n"
    log_fhand.write(msg)
    log_fhand.flush()
    msg = "#STEP 0: converting BAM files to fasta\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    with open(input_sequence) as input_fhand:
        for line in input_fhand:
            if line:
                line = line.rstrip().split()
                group = line[0]
                sub = line[1]
                kind = line[2]
                lowerbound = int(line[3])
                upperbound = int(line[4])
                files = line[-1].split(",")
                checked_files = []
                for file in files:
                    if sequence_kind(file) == "bam":
                        if parser.exclude:
                            msg = "Excluding reads from bed file"
                            log_fhand.write(msg+"\n")
                            results = filter_bam_file(Path(file), Path(parser.output_dir), parser.exclude, parser.num_threads)
                            log = check_run(results)
                            log_fhand.write(log+"\n")
                            log_fhand.flush()
                            print(log)
                            file = results["out_fpath"]
                        results = convert_bam_to_fasta(Path(file), Path(parser.output_dir), parser.num_threads)
                        log = check_run(results)
                        log_fhand.write(log+"\n")
                        log_fhand.flush()
                        print(log)
                        checked_files.append(results["out_fpath"])
                    else:
                        checked_files.append(file)
                data = {"sub": sub, "kind": kind,
                        "lowerbound": lowerbound,
                        "upperbound": upperbound,
                        "files": checked_files}
                if group in inputs:
                    inputs[group].append(data)
                else:
                    inputs[group] = [data]
    return {"inputs": inputs,
            "output": Path(parser.output_dir),
            "threads": parser.num_threads,
            "ram_usage": parser.ram_usage,
            "kmer_size": parser.kmer_size,
            "merge_universe": parser.merge_universe,
            "log": log_fhand}


def main():
    arguments = get_arguments()
    debug = True
    if not debug:
        sys.tracebacklimit = 0
    run_pipeline(arguments)
    arguments["log"].close()


if __name__ == "__main__":
    main()
