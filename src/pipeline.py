from pathlib import Path

from src.kmc import count_kmers, create_input_file, create_kmer_histogram
from src.kmerdata import (load_kmer_table, merge_hetkmers_native, union_universe_size,
                          shannon_estimators_from_counts, kolmogorov_from_table)
from src.kolmogorov import kolmogorov_from_expression_file
from src.utils import check_run, sequence_kind, log_and_print
from src.expression import calculate_sample_estimators as expression_diversity


def build_databases(arguments):
    #STEP 1: build a kmc database per group/sub/dataset, plus one merged
    #database per sub with more than one dataset. Expression-kind datasets
    #are set aside (no kmer counting) for compute_expression_estimators.
    database = {}
    expression = {}
    log_fhand = arguments["log"]
    log_and_print(log_fhand, "#STEP 1: creating databases\n")
    for group, datasets in arguments["inputs"].items():
        database[group] = {}
        occurrences = []
        for dataset in datasets:
            if dataset["kind"] == "expression":
                if group not in expression:
                    expression[group] = {}
                if dataset["sub"] not in expression[group]:
                    expression[group] = {dataset["sub"]: dataset["files"]}
                else:
                    expression[group][dataset["sub"]] += dataset["files"]
                continue
            if dataset["sub"] not in database[group]:
                database[group][dataset["sub"]] = {}
            kinds = [sequence_kind(input_file) for input_file in dataset["files"]]
            if len(set(kinds)) != 1:
                msg = "ERROR: mixed format types found for files {}".format(",".join(dataset["files"]))
                log_and_print(log_fhand, msg)
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
            log_and_print(log_fhand, check_run(results))
            database[group][dataset["sub"]][name] = {"file": results["out_fpath"], "kind": dataset["kind"],
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
            if len(set(kinds)) != 1:
                msg = "ERROR: mixed format types found for files {}".format(",".join(files["files"]))
                log_and_print(log_fhand, msg)
                raise RuntimeError(msg)
            kind = kinds[0]
            input_file_path = create_input_file(files["files"], name, arguments["output"])
            results = count_kmers(input_file_path, name, arguments["output"],
                                    kind, kmer_size=arguments["kmer_size"],
                                    threads=arguments["threads"], max_ram=arguments["ram_usage"])
            log_and_print(log_fhand, check_run(results))
            database[group][sub][name] = {"file": results["out_fpath"], "kind": files["kind"],
                                                    "merged": True}
    return database, expression


def build_histograms(database, log_fhand):
    #STEP 2: kmer histograms, one per database (informational, not consumed
    #by later stages).
    histograms = {}
    log_and_print(log_fhand, "#STEP 2: histograms\n")
    for group, subs in database.items():
        histograms[group] = {}
        for sub, data in subs.items():
            histograms[group][sub] = {}
            for name, values in data.items():
                results = create_kmer_histogram(values["file"], name)
                log_and_print(log_fhand, check_run(results))
                histograms[group][sub][name] = {"file": results["out_fpath"], "kind": values["kind"]}
    return histograms


def dump_counts(database, threads, kmer_size, log_fhand):
    #STEP 3: load raw kmer counts from each database straight into memory
    #(streamed from kmc_tools dump - no .dump file is ever written to disk).
    kmer_tables = {}
    log_and_print(log_fhand, "#STEP 3: creating count dumps\n")
    for group, subs in database.items():
        kmer_tables[group] = {}
        for sub, data in subs.items():
            kmer_tables[group][sub] = {}
            for name, values in data.items():
                #NOTE: reads lowerbound/upperbound from `data` (the per-sub dict),
                #not `values` (the per-name dict) - a known, intentionally
                #unfixed pre-existing bug (see CHANGELOG). Preserved as-is.
                table, results = load_kmer_table(values["file"], kmer_size,
                                            lower_bound=data.get("lowerbound", 1),
                                            upper_bound=data.get("upperbound", 9999999999),
                                            threads=threads, merged=values.get("merged", False))
                table.kind = values["kind"]
                kmer_tables[group][sub][name] = {"table": table, "kind": values["kind"],
                                                    "merged": values.get("merged", False)}
                log_and_print(log_fhand, check_run(results))
    return kmer_tables


def merge_hetkmers(kmer_tables, pooled, log_fhand):
    #STEP 4: for transcriptome data, detect hetkmers (native Hamming-distance-1
    #neighbor search - no smudgeplot) then collapse near-identical kmers
    #together. `pooled` switches between a relaxed component-size cap (pooled
    #multi-individual data) and a strict biallelic + coverage-ratio check
    #(single individual).
    log_and_print(log_fhand, "#STEP 4: native hetkmer detection for transcriptomic data (pooled={})\n".format(pooled))
    for group, subs in kmer_tables.items():
        for sub, data in subs.items():
            for name, values in data.items():
                if values["kind"] != "transcriptome":
                    continue
                merged_table, stats = merge_hetkmers_native(values["table"], pooled=pooled)
                values["table"] = merged_table
                msg = "#SUCCESS: {} components_found={} components_merged={} components_rejected={}".format(
                    name, stats["components_found"], stats["components_merged"], stats["components_rejected"])
                log_and_print(log_fhand, msg)
    return kmer_tables


def compute_universe_sizes(kmer_tables, merge_universe):
    universe_sizes = {}
    if not merge_universe:
        for group, data in kmer_tables.items():
            universe_sizes[group] = {}
            for sub, values in data.items():
                merged = False
                tables = []
                for name, features in values.items():
                    tables.append(features["table"])
                    if features.get("merged", False):
                        merged = True
                        universe_sizes[group][sub] = features["table"].universe_size()
                if not merged:
                    universe_sizes[group][sub] = union_universe_size(tables)
    else:
        for group, data in kmer_tables.items():
            tables_to_combine = []
            for sub, values in data.items():
                merged = False
                tables = []
                for name, features in values.items():
                    tables.append(features["table"])
                    if features.get("merged", False):
                        merged = True
                        tables_to_combine.append(features["table"])
                if not merged:
                    tables_to_combine.extend(tables)
            universe_sizes[group] = union_universe_size(tables_to_combine)
    return universe_sizes


def compute_estimators(kmer_tables, universe_sizes, merge_universe, log_fhand):
    #Shannon diversity + kolmogorov estimators for genomic/transcriptome samples,
    #computed directly from the in-memory KmerTable - no dump/binary files
    #touch disk anywhere in this step. Each sample gets both a regular
    #(count-based) and a presence/absence pass, reported side by side.
    results = {}
    for group, data in kmer_tables.items():
        for sub, values in data.items():
            for name, features in values.items():
                table = features["table"]
                if not merge_universe:
                    universe_size = universe_sizes[group][sub]
                else:
                    universe_size = universe_sizes[group]
                regular = shannon_estimators_from_counts(table.counts, universe_size, presence=False)
                presence = shannon_estimators_from_counts(table.counts, universe_size, presence=True)
                entry = results.setdefault(group, {}).setdefault(sub, {}).setdefault(name, {})
                entry.update(regular)
                entry.update({key + "_presence": value for key, value in presence.items()})
                entry["kolmogorov"] = kolmogorov_from_table(table, universe_size, presence=False)
                entry["kolmogorov_presence"] = kolmogorov_from_table(table, universe_size, presence=True)
                entry.setdefault("universe_size", universe_size)
                entry.setdefault("sub", sub)
                entry.setdefault("name", name)
                entry.setdefault("kind", features["kind"])
                log_and_print(log_fhand, "#SUCCESS: computed estimators for {}".format(name))
    return results


def compute_expression_estimators(expression, results, log_fhand):
    #Shannon diversity + kolmogorov estimators for expression samples, written
    #into the same `results` dict/schema used for genomic/transcriptome samples.
    #Each sample gets both a regular and a presence/absence pass.
    for group, subs in expression.items():
        results.setdefault(group, {})
        for sub, files in subs.items():
            results[group].setdefault(sub, {})
            for count, file in enumerate(files, start=1):
                rep = "{}_{}{}".format(group, sub, count)
                values, universe_size = expression_diversity(Path(file), {}, "TPM", [], binary=False)
                values_presence, _ = expression_diversity(Path(file), {}, "TPM", [], binary=True)
                results[group][sub][rep] = {
                    "kind": "expression", "universe_size": universe_size,
                    "diversity_log2": values["diversity_log2"], "specifity_log2": values["specifity_log2"],
                    "diversity_log10": values["diversity_log10"], "specifity_log10": values["specifity_log10"],
                    "diversity_log2_presence": values_presence["diversity_log2"],
                    "specifity_log2_presence": values_presence["specifity_log2"],
                    "diversity_log10_presence": values_presence["diversity_log10"],
                    "specifity_log10_presence": values_presence["specifity_log10"],
                    "file": Path(file),
                }
                results[group][sub][rep]["kolmogorov"] = kolmogorov_from_expression_file(
                    Path(file), "TPM", [], presence=False)
                results[group][sub][rep]["kolmogorov_presence"] = kolmogorov_from_expression_file(
                    Path(file), "TPM", [], presence=True)
                log_and_print(log_fhand, "#SUCCESS: computed estimators for {}".format(rep))
    return results


def get_files_used(output_dir, prefix):
    filename = Path(output_dir/ "{}.files".format(prefix))
    with open(filename) as fhand:
        return [path.strip() for path in fhand if path]


def write_outputs(results, output_dir):
    with open(output_dir / "file_manifiest.tsv", "w") as manifest_fhand:
        manifest_fhand.write("Group\tSubgroup\tRep\tKind\tFile\n")
        manifest_fhand.flush()
        with open(output_dir / "results.tsv", "w") as out_fhand:
            out_fhand.write("Group\tSubgroup\tRep\tKind\tSubgroup_Universe_Size\t"
                             "Diversity_log2\tSpecifity_log2\tDiversity_log10\tSpecifity_log10\tKolmogorov\t"
                             "Diversity_log2_presence\tSpecifity_log2_presence\tDiversity_log10_presence\t"
                             "Specifity_log10_presence\tKolmogorov_presence\n")
            for group, subs in results.items():
                for sub, reps in subs.items():
                    for rep, features in reps.items():
                        line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                        line = line.format(group, sub, rep, features["kind"],
                                            features["universe_size"], features["diversity_log2"],
                                            features["specifity_log2"], features["diversity_log10"],
                                            features["specifity_log10"], features["kolmogorov"],
                                            features["diversity_log2_presence"], features["specifity_log2_presence"],
                                            features["diversity_log10_presence"], features["specifity_log10_presence"],
                                            features["kolmogorov_presence"])
                        out_fhand.write(line)
                        out_fhand.flush()
                        if features["kind"] == "expression":
                            files_used = [str(features["file"])]
                        else:
                            files_used = get_files_used(output_dir, rep)
                        for file in files_used:
                            fileline = "{}\t{}\t{}\t{}\t{}\n"
                            fileline = fileline.format(group, sub, rep, features["kind"], file)
                            manifest_fhand.write(fileline)
                            manifest_fhand.flush()


def run_pipeline(arguments):
    log_fhand = arguments["log"]
    database, expression = build_databases(arguments)
    histograms = build_histograms(database, log_fhand)
    kmer_tables = dump_counts(database, arguments["threads"], arguments["kmer_size"], log_fhand)
    kmer_tables = merge_hetkmers(kmer_tables, arguments["pooled"], log_fhand)
    universe_sizes = compute_universe_sizes(kmer_tables, arguments["merge_universe"])
    results = compute_estimators(kmer_tables, universe_sizes, arguments["merge_universe"], log_fhand)
    results = compute_expression_estimators(expression, results, log_fhand)
    write_outputs(results, arguments["output"])
    steps = {"database": database, "kmer_tables": kmer_tables,
             "histograms": histograms, "expression": expression}
    return steps, results
