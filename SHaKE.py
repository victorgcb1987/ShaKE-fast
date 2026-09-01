#!/usr/bin/env python


import argparse
import sys

from datetime import datetime
from pathlib import Path
from subprocess import run

from src.kmc import (count_kmers, create_input_file, create_kmer_histogram, 
                     calculate_cutoffs, dump_kmer_counts, calculate_hetkmers)
from src.kolmogorov import calculate_kolmogorov_estimator
from src.utils import check_run, sequence_kind, UnionFind, get_universe_size, convert_bam_to_fasta, filter_bam_file
from src.kmer import calculate_sample_shannon_estimators
from src.expression import calculate_sample_estimators as expression_diversity


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
    help_presence = "(optional) change diversity calculation to presence/absence"
    parser.add_argument("--presence", "-p",
                        help=help_presence, default=False,
                        action="store_true")

    
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
            "presence": parser.presence,
            "log": log_fhand}


def get_files_used(output_dir, prefix):
    filename = Path(output_dir/ "{}.files".format(prefix))
    with open(filename) as fhand:
        return [path.strip() for path in fhand if path]


def main():
    arguments = get_arguments()
    steps = {"database": {}, "count_dumps": {}, "histograms": {}, "expression": {}}
    debug = True
    #Kmer counting
    log_fhand = arguments["log"]
    if not debug:
        sys.tracebacklimit = 0
    arguments = get_arguments()
    
    msg = "#STEP 1: creating databases\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    for group, datasets in arguments["inputs"].items():
        if group not in steps["database"]:
            steps["database"][group] = {}
        occurrences = []
        for dataset in datasets:
            if dataset["kind"] == "expression":
                if group not in steps["expression"]:
                    steps["expression"][group] = {}
                if dataset["sub"] not in steps["expression"][group]:
                    steps["expression"][group] = {dataset["sub"]: dataset["files"]}
                else:
                    steps["expression"][group][dataset["sub"]] += dataset["files"]
                continue
            if dataset["sub"] not in steps["database"][group]:
                steps["database"][group][dataset["sub"]] = {}
            kinds = [sequence_kind(input_file) for input_file in dataset["files"]]
            if len(set(kinds)) !=1:
                msg = "ERROR: mixed format types found for files {}".format(",".join(dataset["files"]))
                log_fhand.write(msg+"\n")
                log_fhand.flush()
                raise RuntimeError(msg)
            kind = kinds[0]
            name = group+"_"+dataset["sub"]
            occurrences.append(name)
            count = str(occurrences.count(name))
            name += count
            input_file_path = create_input_file(dataset["files"], name, arguments["output"])
            results = count_kmers(input_file_path, name, arguments["output"], 
                                    kind, kmer_size=arguments["kmer_size"], 
                                    threads=arguments["threads"], max_ram=arguments["ram_usage"], 
                                    min_occurrence=dataset["lowerbound"], 
                                    max_occurrence=dataset["upperbound"])
            log = check_run(results)
            log_fhand.write(log+"\n")
            log_fhand.flush()
            print(log)
            steps["database"][group][dataset["sub"]][name] = {"file": results["out_fpath"], "kind": dataset["kind"],
                                                                "lowerbound": dataset["lowerbound"], 
                                                                "upperbound": dataset["upperbound"]}
        merges = {}
        for dataset in datasets:
            if dataset["kind"] == "expression":
                continue
            files = [input_file for input_file in dataset["files"]]
            if dataset["sub"] in merges:
                merges[dataset["sub"]]["files"] += files
                merges[dataset["sub"]]["num_datasets"] += 1
            else:
                merges[dataset["sub"]] = {"files": files, 
                                            "num_datasets": 1,
                                            "kind": dataset["kind"]}
        for sub, files in merges.items():
            if files["num_datasets"] == 1:
                continue
            name = group+"_"+sub+"_"+"merged"
            kinds = [sequence_kind(input_file) for input_file in files["files"]]
            if len(set(kinds)) !=1:
                msg = "ERROR: mixed format types found for files {}".format(",".join(dataset["files"]))
                log_fhand.write(msg+"\n")
                log_fhand.flush()
                raise RuntimeError(msg)
            kind = kinds[0]
            input_file_path = create_input_file(files["files"], name, arguments["output"])
            results = count_kmers(input_file_path, name, arguments["output"], 
                                    kind, kmer_size=arguments["kmer_size"], 
                                    threads=arguments["threads"], max_ram=arguments["ram_usage"])
            log = check_run(results)
            log_fhand.write(log+"\n")
            log_fhand.flush()
            print(log)
            steps["database"][group][sub][name] = {"file": results["out_fpath"], "kind": files["kind"], 
                                                    "merged": True}
    msg = "#STEP 2: histograms\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    for group, sub in steps["database"].items():
        if group not in steps["histograms"]:
            steps["histograms"][group] = {}
            for sub, data in sub.items():
                if sub not in steps["histograms"][group]:
                    steps["histograms"][group][sub] = {}
                    for name, values in data.items():
                        results = create_kmer_histogram(values["file"], name)
                        log = check_run(results)
                        log_fhand.write(log+"\n")
                        log_fhand.flush()
                        print(log)
                        steps["histograms"][group][sub][name] = {"file": results["out_fpath"], "kind": values["kind"]}
    msg = "#STEP 3: creating count dumps\n"
    log_fhand.write(msg)
    log_fhand.flush()
    print(msg)
    for group, sub in steps["database"].items():
        if group not in steps["count_dumps"]:
            steps["count_dumps"][group] = {}
            for sub, data in sub.items():
                if sub not in steps["count_dumps"][group]:
                    steps["count_dumps"][group][sub] = {}
                    for name, values in data.items():
                        results = dump_kmer_counts(values["file"], name,
                                                    lower_bound=data.get("lowerbound", 1),
                                                    upper_bound=data.get("upperbound", 9999999999),
                                                    threads=arguments["threads"], pipe=True)
                        steps["count_dumps"][group][sub][name] = {"file": results["out_fpath"], "kind": values["kind"],
                                                                    "merged": values.get("merged", False),
                                                                    }
                        log = check_run(results)
                        log_fhand.write(log+"\n")
                        log_fhand.flush()
                        print(log)
    msg = "#STEP 4: calculating hetkmers for transcriptomic data\n"
    print(msg)
    log_fhand.write(msg)
    log_fhand.flush()
    for group, sub in steps["count_dumps"].items():
        for sub, data in sub.items():
            for name, values in data.items():
                if values["kind"] == "transcriptome":
                    results = calculate_hetkmers(values["file"], arguments["output"])
                    log = check_run(results)
                    print(log)
                    log_fhand.write(log+"\n")
                    log_fhand.flush()
                    merged_by_hetkmer = Path(arguments["output"] / "{}_grouped_by_hetkmers.dump".format(values["file"].name.replace(".dump", "")))
                    if merged_by_hetkmer.exists():
                        msg = "#ALREADY_DONE: merged {} into {}\n".format(str(values["file"]), str(merged_by_hetkmer))
                        print(msg)
                        log_fhand.write(msg)
                        log_fhand.flush()
                        values["file"] = Path(merged_by_hetkmer)
                        continue

                    with open(results["out_fpath"]) as fhand:
                        hetkmers = [(line.rstrip().split()[1], line.split()[0]) for line in fhand]
                    with open(values["file"]) as fhand:
                        values_ = {line.split()[0]: int(line.rstrip().split()[1]) for line in fhand}
                    #union-find algorithm
                    unique_elements = set(values_.keys())
                    for hetkmer in hetkmers:
                        unique_elements.update(hetkmer)
                    uf = UnionFind(unique_elements)
                    for hetkmer in hetkmers:
                        uf.join(hetkmer[0], hetkmer[1])
                    groups = {}
                    for seq, value in values_.items():
                        representative = uf.find(seq)
                        if representative in groups:
                            groups[representative] += value
                        else:
                            groups[representative] = value
                    with open(merged_by_hetkmer, "w") as out_fhand:
                        for seq, value in groups.items():
                            out_fhand.write(f'{seq}\t{value}\n')
                    values["file"] = Path(merged_by_hetkmer)
        
    universes_sizes = {}
    if not arguments["merge_universe"]:
        for group, data in steps["count_dumps"].items():
            if group not in universes_sizes:
                universes_sizes[group] = {}
            for sub, values in data.items():
                merged = False
                files = []
                for name, features in values.items():
                    files.append(str(features["file"]))
                    if features.get("merged", False):
                        merged = True
                        universe_size = get_universe_size([str(features["file"])])
                        universes_sizes[group][sub] = universe_size
                if not merged:
                    universe_size = get_universe_size(files)
                    universes_sizes[group][sub] = universe_size
    else:
        for group, data in steps["count_dumps"].items():
            if group not in universes_sizes:
                universes_sizes[group] = []
                files_to_combine = []
            for sub, values in data.items():
                merged = False
                files = []
                for name, features in values.items():
                    files.append(str(features["file"]))
                    if features.get("merged", False):
                        merged = True
                        files_to_combine.append(str(features["file"]))
                if not merged:
                    for file in files:
                        files_to_combine.append(file)
            universe_size = get_universe_size(files_to_combine)
            universes_sizes[group] = universe_size            
        
    steps["estimators"] = {}
    results = {}
    print(steps["count_dumps"])
    for group, data in steps["count_dumps"].items():
        if group not in  steps["estimators"]:
            steps["estimators"][group] = {}
            for sub, values in data.items():
                for name, features in values.items():
                    if not arguments["merge_universe"]:
                        universe_size = universes_sizes[group][sub]
                    else:
                        universe_size = universes_sizes[group]
                    print(features["file"])
                    calculate_sample_shannon_estimators(features["file"], universe_size, results, group=group, 
                                                        sub=sub, name=name, kind=features["kind"], file=features["file"], 
                                                        pipe=True, binary=arguments["presence"])
                    calculate_kolmogorov_estimator(features["file"], universe_size, results, group=group, 
                                                    sub=sub, name=name, kind=features["kind"], units="TPM",
                                                    presence=arguments["presence"])
                        
    with open(arguments["output"] / "file_manifiest.tsv", "w") as manifiest_fhand:
        manifiest_fhand.write("Group\tSubgroup\tRep\tKind\"File\n")
        manifiest_fhand.flush()
        with open(arguments["output"] / "results.tsv", "w") as out_fhand:
            out_fhand.write("Group\tSubgroup\tRep\tKind\tSubgroup_Universe_Size\tDiversity_log2\tSpecifity_log2\tDiversity_log10\tSpecifity_log10\tKolmogorov\n")
            for group, values in results.items():
                for sub, reps in values.items():
                    for rep, features in reps.items():
                        line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                        line = line.format(group, sub, rep, features["kind"], 
                                            features["universe_size"], features["diversity_log2"],
                                            features["specifity_log2"], features["diversity_log10"],
                                            features["specifity_log10"], features["kolmogorov"])
                        out_fhand.write(line)
                        out_fhand.flush()
                        files_used = get_files_used(arguments["output"], rep)
                        for file in files_used:
                            fileline = "{}\t{}\t{}\t{}\t{}\n"
                            fileline = fileline.format(group, sub, rep, features["kind"], file)
                            manifiest_fhand.write(fileline)
                            manifiest_fhand.flush()    

         
            expression_estimators = {}
            for group,  subs in steps["expression"].items():
                for sub, files in subs.items():
                    count = 0
                    for file in files:
                        count += 1
                        rep = "{}_{}{}".format(group, sub,count)
                        values, universe_size = expression_diversity(Path(file), expression_estimators, "TPM", [], 
                                                                        binary=arguments["presence"])
                        calculate_kolmogorov_estimator(features["file"], universe_size, results, group=group, 
                                                    sub=sub, name=name, kind="expression", units="TPM")
                        
                        line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                        line = line.format(group, sub, rep, "expression", universe_size,
                                           features["diversity_log2"], features["specifity_log2"], 
                                           features["diversity_log10"], features["specifity_log10"], 
                                           features["kolmogorov"])
                        out_fhand.write(line)
                        out_fhand.flush()
                        fileline = "{}\t{}\t{}\t{}\t{}\n"
                        fileline = fileline.format(group, sub, rep, "expression", file)
                        manifiest_fhand.write(fileline)
                        manifiest_fhand.flush()
    log_fhand.close()


if __name__ == "__main__":
    main()
