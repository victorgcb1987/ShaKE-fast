#!/usr/bin/env python


import argparse
import sys

from pathlib import Path

from src.kmc import calculate_hetkmers

def parse_arguments():
    desc = "get hetkmers between one or more kmc count dumps"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = "(Required) directory with dumps. dumps should have .dump filename extension"
    parser.add_argument("--input_file", "-i", type=str,
                        help=help_input, required=True)
    
    help_output = "(Required) output file"
    parser.add_argument("--output_dir", "-o", type=str,
                        help=help_output, required=True)
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    input_fdir = Path(parser.input_file)
    output_fpath = Path(parser.output_dir)
    return {"input_folder": input_fdir,
            "out_fpath": output_fpath}


def main():
    arguments = get_arguments()
    hetkmer_files = [calculate_hetkmers(dump_file, arguments["out_fpath"])["out_fpath"] for dump_file in arguments["input_folder"].glob("*.dump")]

if __name__ == "__main__":
    main()