#!/usr/bin/env python


import argparse
import sys

from datetime import datetime
from pathlib import Path

from src.kmc import (count_kmers, create_input_file, create_kmer_histogram, 
                     calculate_cutoffs, dump_kmer_counts)
from src.utils import check_run, sequence_kind


def parse_arguments():
    desc = "Calculate kmers occurrences from fastq/fasta files"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = '''(Required) TSV file containing accession, kind 
                    and files separated by commas, for example
                    ACC1    file1,file2
                    ACC2    file3
                    If there is more than one file for accession
                    kmers while be combined by union-sum'''
    parser.add_argument("--input_file", "-i", type=str,
                        help=help_input, required=True)
    
    help_output = "(Required) output dir"
    parser.add_argument("--output_dir", "-o", type=str,
                        help=help_output, required=True)
    
    help_kmer_length = "(Required) kmer length"
    parser.add_argument("--kmer_length", "-k", type=int,
                        help=help_kmer_length, required=True)
    
    help_lower_bound = "(Optional) lower occurrence in kmers"
    parser.add_argument("--lower_bound", "-l", type=int,
                        help=help_lower_bound, required=False,
                        default=0)
    
    help_upper_bound = "(Optional) upper occurrence in kmers"
    parser.add_argument("--upper_bound", "-u", type=int,
                        help=help_upper_bound, required=False,
                        default=10000000000)
    help_calculate_bounds = "(optional) calculate bounds"
    parser.add_argument("--calculate_bounds", "-c", action="store_true",
                        help=help_calculate_bounds)

    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    inputs = {}
    input_sequence = Path(parser.input_file)
    with open(input_sequence) as input_fhand:
        universe_files = []
        for line in input_fhand:
            if line:
                files = line.split()[1].split(",")
                for _file in files:
                    if _file not in universe_files:
                        universe_files.append(_file)
        inputs["universe"] = universe_files        
    return {"inputs": inputs,
            "output": Path(parser.output_dir),
            "kmer_length": parser.kmer_length,
            "lower_bound": parser.lower_bound,
            "upper_bound": parser.upper_bound,
            "calculate_bounds": parser.calculate_bounds}

def main():
    arguments = get_arguments()
    lower_bound = arguments["lower_bound"]
    upper_bound = arguments["upper_bound"]
    logdate = "Kmer_counting_"+datetime.now().strftime("%d_%m_%Y-%H_%M_%S") + ".log"
    log_fname = arguments["output"] / logdate
    if not arguments["output"].exists():
        arguments["output"].mkdir(exist_ok=True, parents=True)
        
    with open(log_fname, "w") as log_fhand:
        log_fhand.write("#Command Used: "+ "".join(sys.argv)+"\n")
        for name, input_files in arguments["inputs"].items():
            kinds = [sequence_kind(input_file) for input_file in input_files]
            if len(set(kinds)) != 1:
                raise ValueError("mixed format types for {} found".format(name))
            kind = kinds[0]
            input_file_path = create_input_file(input_files, name, 
                                                arguments["output"])
            results = count_kmers(input_file_path, name, 
                                  arguments["output"], kind, kmer_size=21, 
                                  threads=40, max_ram=64, min_occurrence=lower_bound,
                                  max_occurrence=upper_bound)
            log_fhand.write(check_run(results)+"\n")
            if arguments["calculate_bounds"]:
                lower_bound, upper_bound = calculate_cutoffs(results["out_fpath"])
            results = dump_kmer_counts(Path(arguments["output"]), name, threads=40, 
                                       lower_bound=lower_bound, upper_bound=upper_bound)
            log_fhand.write(check_run(results)+"\n")
  
if __name__ == "__main__":
    main()

