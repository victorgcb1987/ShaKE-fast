#!/usr/bin/env python

import argparse
import sys

from pathlib import Path

from src.kmer import calculate_sample_estimators


def parse_arguments():
    desc = "Create a table with "
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = '''(Required) Directory with meryl counts.
                    Subolders have the suffix .count'''
    parser.add_argument("--input_dir", "-i", type=str,
                        help=help_input, required=True)
    
    help_output = "(Required) output dir"
    parser.add_argument("--output", "-o", type=str,
                        help=help_output, required=True)
    help_output = "(Optional) write kmer_counts_into a file"
    parser.add_argument("--write_kmers", "-w", action="store_true",
                        default=False)
    
    help_output = "(Required) Universe size"
    parser.add_argument("--universe_size", "-N", required=True,
                        type=int)
    help_binary = "(Optioanl) change expression values to 1 or 0. False by default"
    parser.add_argument("-b", "--binary", action="store_true",
                        default=False, help=help_binary)

    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    inputs = [dir for dir in Path(parser.input_dir).glob("*.dump")]
    return {"inputs": inputs,
            "output": Path(parser.output),
            "write_kmers": parser.write_kmers,
            "N": parser.universe_size,
            "binary": parser.binary}

def main():
    arguments = get_arguments()
    output_dir = arguments["output"]
    if not output_dir.exists():
        output_dir.mkdir(exist_ok=True, parents=True)
    estimators = {}
    N = arguments["N"]
    for filepath in arguments["inputs"]:
         calculate_sample_estimators(filepath, N, estimators, binary=arguments["binary"])
    with open(output_dir/"sample_estimators.tsv", "w") as out_fhand:
        out_fhand.write("Sample\tDiversity\tSpecificity\n")
        for sample, values in estimators.items():
            out_fhand.write("{}\t{}\t{}\n".format(sample, values["diversity"], values["specifity"]))

if __name__ == "__main__":
    main()