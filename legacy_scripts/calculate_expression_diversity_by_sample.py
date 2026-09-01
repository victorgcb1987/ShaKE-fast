#!/usr/bin/env python

import argparse
import sys

from pathlib import Path

from src.expression import calculate_sample_estimators


def parse_arguments():
    desc = "Calculate expression "
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = '''(Required) Directory with transcript/gene abundancy counts.
                 count files should have the prefix .abund.tsv'''
    parser.add_argument("--input_dir", "-i", type=str,
                        help=help_input, required=True)
    help_units = "Measurement units. TPM by default"
    parser.add_argument("--units", "-u", type=str,
                        help=help_units, default="TPM")
    help_output = "(Required) output dir"
    parser.add_argument("--output", "-o", type=str,
                        help=help_output, required=True)
    help_exclude = "(Optional) assembly sequences to exclude"
    parser.add_argument('-e', '--exclude', nargs='+', default=[],
                        help=help_exclude)
    help_binary = "(Optioanl) change expression values to 1 or 0. False by default"
    parser.add_argument("-b", "--binary", action="store_true",
                        default=False, help=help_binary)
    

    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()



def get_arguments():
    parser = parse_arguments()
    inputs = [dir for dir in Path(parser.input_dir).glob("*.abund.tsv")]
    return {"inputs": inputs,
            "output": Path(parser.output),
            "units": parser.units,
            "exclude": parser.exclude,
            "binary": parser.binary}


def main():
    arguments = get_arguments()
    output_dir = arguments["output"]
    if not output_dir.exists():
        output_dir.mkdir(exist_ok=True, parents=True)
    estimators = {}
    for filepath in arguments["inputs"]:
         calculate_sample_estimators(filepath, estimators, arguments["units"], arguments["exclude"], 
                                     binary=arguments["binary"])
    with open(output_dir/"sample_estimators.tsv", "w") as out_fhand:
        out_fhand.write("Sample\tDiversity\tSpecificity\n")
        for sample, values in estimators.items():
            out_fhand.write("{}\t{}\t{}\n".format(sample, values["diversity"], values["specifity"]))


if __name__ == "__main__":
    main()