#!/usr/bin/env python

import argparse
import multiprocess as mp
import sys
from math import log10 as log


from itertools import islice
from pathlib import Path
from shutil import rmtree as remove_dir
from src.utils import merge_temporary_files, get_kmer_value


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
    
    help_output = "(Required) Dump file with universe of kmers"
    parser.add_argument("--universe", "-u", required=True,
                        type=str)
    help_threads = "(Optional) number of threads. 6 by default"
    parser.add_argument("--threads", "-t", default=6,
                       type=int)

    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    inputs = sorted([dir for dir in Path(parser.input_dir).glob("*.dump")])
    return {"inputs": inputs,
            "output": Path(parser.output),
            "universe": Path(parser.universe),
            "N": len(inputs),
            "threads": parser.threads}


def chunk_processing(generator, chunk_size):
    iterator = iter(generator)
    while True:
        chunk = list(islice(iterator, chunk_size))
        if not chunk:
            break
        yield chunk


def calculate_kmer_estimators(kmer, filepaths, universe_size):
     raw_values = [get_kmer_value(filepath, kmer) for filepath in filepaths]
     N = sum(raw_values)
     values = [(float(value)/N) * log(float(value)/N) if value > 0 else 0 for value in raw_values]
     diversity_value =  -sum(value for value in values if value != 0)
     specifity = log(universe_size) - diversity_value
     return kmer, diversity_value, specifity 


def main():
    arguments = get_arguments()
    tmp_dir = arguments["output"] / "tmp"
    tmp_dir.mkdir(parents=False, exist_ok=True)
    with open(tmp_dir / "Header.estim", "w") as header_fhand:
        header_fhand.write("Kmer\tDiversity\tSpecifity\n")
    k = 0
    with open(arguments["universe"]) as universe_fhand:
        kmer_list = (line.split()[0] for line in universe_fhand)
        def pack_args(kmer):
            return calculate_kmer_estimators(kmer, arguments["inputs"], len(arguments["inputs"]))
        for chunk in chunk_processing(kmer_list, arguments["threads"]):
            with mp.Pool(processes=arguments["threads"]) as pool:
                results = pool.map(pack_args, chunk)
                with open(tmp_dir/"K{}.estim".format(k), "w") as chunk_fhand:
                    with open(tmp_dir/"K{}.index".format(k), "w") as index_fhand:
                        for result in results:
                            k += 1
                            chunk_fhand.write("K{}\t{}\t{}\n".format(k, result[1], result[2]))
                            index_fhand.write("K{}\t{}\n".format(k, result[0]))

        merge_temporary_files(tmp_dir, arguments["output"], ".estim")
        merge_temporary_files(tmp_dir, arguments["output"], ".index")
        remove_dir(tmp_dir)

if __name__ == "__main__":
    main()