from pathlib import Path

from src.kmc import (count_kmers, create_input_file, create_kmer_histogram,
                     remove_kmc_database)
import numpy as np

from itertools import groupby

from src.kmerdata import (load_kmer_table, merge_hetkmers_native, union_universe_size,
                          shannon_estimators_from_counts, kolmogorov_estimators_from_counts,
                          shannon_diversity_log10, shannon_from_diversity_log10,
                          evenness_from_diversity_log10,
                          StreamingKolmogorov, kolmogorov_estimators_from_streams)
from src.kolmogorov import kolmogorov_from_expression_file
from src.utils import check_run, sequence_kind, log_and_print
from src.expression import calculate_sample_estimators as expression_diversity


def plan_jobs(arguments):
    #Works out, without running anything, every kmc database that has to be
    #built: one per (non-expression) dataset, plus one merged database per sub
    #with more than one dataset. Returns (jobs, expression); `jobs` is in the
    #order they must be built (per group: datasets first, then merged subs).
    #Expression-kind datasets are set aside (no kmer counting) for
    #compute_expression_estimators.
    log_fhand = arguments["log"]
    jobs = []
    expression = {}
    for group, datasets in arguments["inputs"].items():
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
            kinds = [sequence_kind(input_file) for input_file in dataset["files"]]
            if len(set(kinds)) != 1:
                msg = "ERROR: mixed format types found for files {}".format(",".join(dataset["files"]))
                log_and_print(log_fhand, msg)
                raise RuntimeError(msg)
            name = group+"_"+dataset["sub"]
            occurrences.append(name)
            name += str(occurrences.count(name))
            jobs.append({"group": group, "sub": dataset["sub"], "name": name,
                         "files": dataset["files"], "seq_kind": kinds[0],
                         "kind": dataset["kind"], "merged": False,
                         "lowerbound": dataset["lowerbound"],
                         "upperbound": dataset["upperbound"]})
        merges = {}
        for dataset in datasets:
            if dataset["kind"] == "expression":
                continue
            files = [input_file for input_file in dataset["files"]]
            bounds = (dataset["lowerbound"], dataset["upperbound"])
            if dataset["sub"] in merges:
                merges[dataset["sub"]]["files"] += files
                merges[dataset["sub"]]["num_datasets"] += 1
                merges[dataset["sub"]]["bounds"].add(bounds)
            else:
                merges[dataset["sub"]] = {"files": files,
                                            "num_datasets": 1,
                                            "kind": dataset["kind"],
                                            "bounds": {bounds}}
        for sub, files in merges.items():
            if files["num_datasets"] == 1:
                continue
            kinds = [sequence_kind(input_file) for input_file in files["files"]]
            if len(set(kinds)) != 1:
                msg = "ERROR: mixed format types found for files {}".format(",".join(files["files"]))
                log_and_print(log_fhand, msg)
                raise RuntimeError(msg)
            job = {"group": group, "sub": sub, "name": group+"_"+sub+"_"+"merged",
                   "files": files["files"], "seq_kind": kinds[0],
                   "kind": files["kind"], "merged": True}
            #The merged database gets the datasets' cutoffs only if they all
            #agree; otherwise there is no single cutoff to apply to it.
            if len(files["bounds"]) == 1:
                job["lowerbound"], job["upperbound"] = next(iter(files["bounds"]))
            else:
                msg = ("#WARNING: datasets of {} {} have different count cutoffs {}; "
                       "the merged database is built without cutoffs".format(
                           group, sub, sorted(files["bounds"])))
                log_and_print(log_fhand, msg)
            jobs.append(job)
    #A sub's universe is its merged database if it has one, else its only
    #dataset (a sub with several datasets always has a merged one) - this
    #flags the jobs that define it.
    with_merged = {(job["group"], job["sub"]) for job in jobs if job["merged"]}
    for job in jobs:
        job["in_universe"] = job["merged"] or (job["group"], job["sub"]) not in with_merged
    return jobs, expression


def count_job(job, arguments):
    #Builds the kmc database of one planned job (see plan_jobs).
    input_file_path = create_input_file(job["files"], job["name"], arguments["output"])
    cutoffs = {}
    if "lowerbound" in job:
        cutoffs = {"min_occurrence": job["lowerbound"], "max_occurrence": job["upperbound"]}
    return count_kmers(input_file_path, job["name"], arguments["output"],
                       job["seq_kind"], kmer_size=arguments["kmer_size"],
                       threads=arguments["threads"], max_ram=arguments["ram_usage"],
                       **cutoffs)


def build_databases(arguments):
    #STEP 1: build a kmc database per group/sub/dataset, plus one merged
    #database per sub with more than one dataset.
    log_fhand = arguments["log"]
    log_and_print(log_fhand, "#STEP 1: creating databases\n")
    jobs, expression = plan_jobs(arguments)
    database = {}
    for group, datasets in arguments["inputs"].items():
        database[group] = {}
        for dataset in datasets:
            if dataset["kind"] != "expression":
                database[group].setdefault(dataset["sub"], {})
    for job in jobs:
        results = count_job(job, arguments)
        log_and_print(log_fhand, check_run(results))
        entry = {"file": results["out_fpath"], "kind": job["kind"]}
        if job["merged"]:
            entry["merged"] = True
        if "lowerbound" in job:
            entry.update({"lowerbound": job["lowerbound"], "upperbound": job["upperbound"]})
        database[job["group"]][job["sub"]][job["name"]] = entry
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
                table, results = load_kmer_table(values["file"], kmer_size,
                                            lower_bound=values.get("lowerbound", 1),
                                            upper_bound=values.get("upperbound", 9999999999),
                                            threads=threads, merged=values.get("merged", False))
                table.kind = values["kind"]
                kmer_tables[group][sub][name] = {"table": table, "kind": values["kind"],
                                                    "merged": values.get("merged", False)}
                log_and_print(log_fhand, check_run(results))
    return kmer_tables


def merge_sample_hetkmers(name, values, pooled, log_fhand):
    #Collapses near-identical kmers of one transcriptome sample in place
    #(values is the {"table", "kind", ...} entry of kmer_tables).
    if values["kind"] != "transcriptome":
        return
    merged_table, stats = merge_hetkmers_native(values["table"], pooled=pooled)
    values["table"] = merged_table
    msg = "#SUCCESS: {} components_found={} components_merged={} components_rejected={}".format(
        name, stats["components_found"], stats["components_merged"], stats["components_rejected"])
    log_and_print(log_fhand, msg)


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
                merge_sample_hetkmers(name, values, pooled, log_fhand)
    return kmer_tables


def process_samples_sequentially(arguments):
    #Low-disk/low-memory alternative to STEPS 1-6: every sample goes through
    #its whole chain (kmc database -> histogram -> in-memory kmer table ->
    #het-kmer merge) on its own and is then reduced to what the results need;
    #its kmc database (.kmc_pre/.kmc_suf, the bulk of the disk usage) and its
    #kmer table are dropped before the next sample starts. Only a database
    #built by this run is deleted - one that was already there is left alone.
    #
    #What a sample keeps: its Shannon diversity (it doesn't depend on the
    #universe size), and a streaming Kolmogorov compressor already fed with
    #its own counts, waiting for the zeros that the universe size will add.
    #The universe size depends on the other samples, so specificity and the
    #Kolmogorov ratio are finished once the last job of the group is done;
    #until then only the set of kmers defining the universe is kept (a single
    #running union per group with --merge_universe, nothing at all otherwise,
    #since a sub's universe is just the size of its merged/only table).
    log_fhand = arguments["log"]
    log_and_print(log_fhand, "#SEQUENTIAL MODE: kmc databases and kmer tables are dropped as soon as each sample is summarized\n")
    jobs, expression = plan_jobs(arguments)
    results = {}
    for group, group_jobs in groupby(jobs, key=lambda job: job["group"]):
        group_jobs = list(group_jobs)
        samples = []
        sub_universe = {}
        group_kmers = np.empty(0, dtype=np.uint64)
        for job in group_jobs:
            name = job["name"]
            log_and_print(log_fhand, "#SAMPLE: {}".format(name))
            count_results = count_job(job, arguments)
            log_and_print(log_fhand, check_run(count_results))
            db_fpath = count_results["out_fpath"]
            log_and_print(log_fhand, check_run(create_kmer_histogram(db_fpath, name)))
            table, dump_results = load_kmer_table(db_fpath, arguments["kmer_size"],
                                                  lower_bound=job.get("lowerbound", 1),
                                                  upper_bound=job.get("upperbound", 9999999999),
                                                  threads=arguments["threads"], merged=job["merged"])
            table.kind = job["kind"]
            log_and_print(log_fhand, check_run(dump_results))
            if dump_results["returncode"] == 0 and count_results["returncode"] == 0:
                removed = remove_kmc_database(db_fpath)
                log_and_print(log_fhand, "#CLEANED: {} ({} files removed)".format(name, len(removed)))
            values = {"table": table, "kind": job["kind"], "merged": job["merged"]}
            merge_sample_hetkmers(name, values, arguments["pooled"], log_fhand)
            table = values["table"]
            if job["in_universe"]:
                if arguments["merge_universe"]:
                    group_kmers = np.union1d(group_kmers, table.kmers)
                else:
                    sub_universe[job["sub"]] = table.universe_size()
            samples.append({"job": job,
                            "diversity_log10": shannon_diversity_log10(table.counts, presence=False),
                            "diversity_log10_presence": shannon_diversity_log10(table.counts, presence=True),
                            "kolmogorov": StreamingKolmogorov(table.counts, presence=False),
                            "kolmogorov_presence": StreamingKolmogorov(table.counts, presence=True),
                            "total_count": int(table.counts.sum()),
                            "n_positive": int((table.counts > 0).sum())})
            del table, values
        group_universe = int(group_kmers.shape[0])
        del group_kmers
        #same sub order as the regular flow: first appearance among the datasets
        for sub in dict.fromkeys(job["sub"] for job in group_jobs):
            for sample in samples:
                job = sample["job"]
                if job["sub"] != sub:
                    continue
                universe_size = group_universe if arguments["merge_universe"] else sub_universe[sub]
                entry = results.setdefault(group, {}).setdefault(sub, {}).setdefault(job["name"], {})
                entry.update(shannon_from_diversity_log10(sample["diversity_log10"], universe_size))
                entry.update({key + "_presence": value for key, value in
                              shannon_from_diversity_log10(sample["diversity_log10_presence"], universe_size).items()})
                entry.update(evenness_from_diversity_log10(sample["diversity_log10"],
                                                           sample["n_positive"], universe_size))
                entry.update(kolmogorov_estimators_from_streams(
                    sample["kolmogorov"], sample["kolmogorov_presence"],
                    sample["total_count"], universe_size))
                entry.setdefault("universe_size", universe_size)
                entry.setdefault("sub", sub)
                entry.setdefault("name", job["name"])
                entry.setdefault("kind", job["kind"])
                log_and_print(log_fhand, "#SUCCESS: computed estimators for {}".format(job["name"]))
    tmp_dir = arguments["output"] / "tmp"
    if tmp_dir.is_dir() and not any(tmp_dir.iterdir()):
        tmp_dir.rmdir()
    return results, expression


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
                entry.update(evenness_from_diversity_log10(regular["diversity_log10"],
                                                           int((table.counts > 0).sum()), universe_size))
                entry.update(kolmogorov_estimators_from_counts(table.counts, universe_size))
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
                results[group][sub][rep].update(evenness_from_diversity_log10(
                    values["diversity_log10"], values["n_positive"], universe_size))
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
                             "Specifity_log10_presence\tKolmogorov_presence\tKolmogorov_norm\t"
                             "Pielou_Evenness\tUniverse_Evenness\n")
            for group, subs in results.items():
                for sub, reps in subs.items():
                    for rep, features in reps.items():
                        line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                        line = line.format(group, sub, rep, features["kind"],
                                            features["universe_size"], features["diversity_log2"],
                                            features["specifity_log2"], features["diversity_log10"],
                                            features["specifity_log10"], features["kolmogorov"],
                                            features["diversity_log2_presence"], features["specifity_log2_presence"],
                                            features["diversity_log10_presence"], features["specifity_log10_presence"],
                                            features["kolmogorov_presence"],
                                            features.get("kolmogorov_norm", "NA"),
                                            features["pielou_evenness"], features["universe_evenness"])
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
    if arguments.get("sequential", False):
        #estimators are computed along the way - no tables are kept
        results, expression = process_samples_sequentially(arguments)
        steps = {"database": {}, "kmer_tables": {}, "histograms": {}, "expression": expression}
    else:
        database, expression = build_databases(arguments)
        histograms = build_histograms(database, log_fhand)
        kmer_tables = dump_counts(database, arguments["threads"], arguments["kmer_size"], log_fhand)
        kmer_tables = merge_hetkmers(kmer_tables, arguments["pooled"], log_fhand)
        universe_sizes = compute_universe_sizes(kmer_tables, arguments["merge_universe"])
        results = compute_estimators(kmer_tables, universe_sizes, arguments["merge_universe"], log_fhand)
        steps = {"database": database, "kmer_tables": kmer_tables,
                 "histograms": histograms, "expression": expression}
    results = compute_expression_estimators(expression, results, log_fhand)
    write_outputs(results, arguments["output"])
    return steps, results
