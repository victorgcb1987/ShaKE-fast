#!/usr/bin/env python


import argparse
import sys

from pathlib import Path

from src.kmc import calculate_hetkmers
from src.utils import  UnionFind


def parse_arguments():
    desc = "get hetkmers between one or more kmc count dumps"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = "(Required) directory with dumps. dumps should have .dump filename extension"
    parser.add_argument("--input_dir", "-i", type=str,
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
    input_fdir = Path(parser.input_dir)
    output_fpath = Path(parser.output_dir)
    return {"input_dir": input_fdir,
            "out_fpath": output_fpath}

def main():
    arguments = get_arguments()
    input_folder = arguments["input_dir"]
    for fpath in input_folder.glob("*.dump"):
        hetkmer_seqs = input_folder / "{}_hetkmers_sequences.tsv".format(fpath.name)
        if not hetkmer_seqs.exists():
            print("hetkmer not found: {}".format(str(hetkmer_seqs)))
        else:
            with open(hetkmer_seqs) as fhand:
                hetkmers = [(line.rstrip().split()[1], line.split()[0]) for line in fhand]
            with open(fpath) as fhand:
                values = {line.split()[0]: int(line.rstrip().split()[1]) for line in fhand}
        #union-find algorithm
        unique_elements = set(values.keys())
        for hetkmer in hetkmers:
            unique_elements.update(hetkmer)
        uf = UnionFind(unique_elements)
        for hetkmer in hetkmers:
            uf.join(hetkmer[0], hetkmer[1])

        groups = {}
        for seq, value in values.items():
            representative = uf.find(seq)
            if representative in groups:
                groups[representative] += value
            else:
                groups[representative] = value
        output_fpath = arguments["out_fpath"] / "{}_grouped_by_hetkmers.dump".format(fpath.name.replace(".dump", ""))
        with open(output_fpath, "w") as out_fhand:
            for seq, value in groups.items():
                out_fhand.write(f'{seq}\t{value}\n')        

if __name__ == "__main__":
    main()

