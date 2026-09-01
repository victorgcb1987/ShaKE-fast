from subprocess import run
from pathlib import Path


def create_input_file(fpaths, name, output_fpath):
    filepats_fpath = output_fpath / "{}.files".format(name)
    with open(filepats_fpath, "w") as filepaths_fhand:
        for fpath in fpaths:
            filepaths_fhand.write(str(fpath)+"\n")
    return filepats_fpath


def count_kmers(input_file, name, output_dir, kind, kmer_size=21, 
                threads=6, max_ram=6, min_occurrence=1, 
                max_occurrence=10000000000):
    temp_dir = output_dir/"tmp"
    if not temp_dir.exists():
        temp_dir.mkdir(parents=True)
    out_db_fpath = output_dir / name
    cmd = "kmc "
    if kind == "fasta":
        cmd += "-fm "
    elif kind == "fastq":
        cmd += "-fq "
    cmd += "-k{} -t{} -m{} -sm -ci{} -cx{} -cs10000000000 @{} {} {} "
    cmd = cmd.format(kmer_size, threads, max_ram, 
                     min_occurrence, max_occurrence, str(input_file),
                    str(out_db_fpath), str(temp_dir))
    out_pre = Path("{}.kmc_pre".format(out_db_fpath))
    out_suf = Path("{}.kmc_suf".format(out_db_fpath))
    if out_pre.exists() and out_suf.exists():
        results = {"command": cmd, "returncode": 99, "name": name,
                   "msg": "output file exists already", "out_fpath": out_db_fpath}
    else:
        run_ = run(cmd, shell=True, capture_output=True)
        results = {"command": cmd, "returncode": run_.returncode, "name": name,
                   "msg": run_.stderr.decode(), "out_fpath": out_db_fpath}
    return results
    

def create_kmer_histogram(input_db_fpath, name):
    out_histogram_fpath = "{}.hist".format(str(input_db_fpath))
    cmd = "kmc_tools transform {} histogram {}".format(input_db_fpath,
                                                       out_histogram_fpath)
    if Path(out_histogram_fpath).exists():
        results = {"command": cmd, "returncode": 99, "name": name,
                   "msg": "output file exists already", 
                   "out_fpath": out_histogram_fpath}
    else:
        run_ = run(cmd, shell=True, capture_output=True)
        results = {"command": cmd, "returncode": run_.returncode, "name": name,
               "msg": run_.stderr.decode(), "out_fpath": out_histogram_fpath} 
    return results


def calculate_cutoffs(histogram_fpath):
    lower_cmd = "smudgeplot.py cutoff {} L".format(str(histogram_fpath))
    run_ = run(lower_cmd, shell=True, capture_output=True)
    lower_bound = run_.stdout.decode()
    upper_cmd = "smudgeplot.py cutoff {} U".format(str(histogram_fpath))
    run_ = run(upper_cmd, shell=True, capture_output=True)
    upper_bound = run_.stdout.decode()
    return lower_bound, upper_bound


def dump_kmer_counts(db_fpath, name, threads=6, lower_bound=1, upper_bound=99999999999, pipe=False):
    if not pipe:
        db_fpath = str(db_fpath / name)
    out_dump_fpath = "{}_L{}_U{}.dump"
    out_dump_fpath = out_dump_fpath.format(db_fpath, lower_bound, 
                                           upper_bound)
    cmd = "kmc_tools -t{} transform {} -ci{} -cx{} dump -s {}"
    cmd = cmd.format(threads, db_fpath, lower_bound, 
                     upper_bound, out_dump_fpath)
    if Path(out_dump_fpath).exists():
        results = {"command": cmd, "returncode": 99, "name": name,
                   "msg": "output file exists already", 
                   "out_fpath": Path(out_dump_fpath)}
    else:
        run_ = run(cmd, shell=True, capture_output=True)
        results = {"command": cmd, "returncode": run_.returncode, "name": name,
                   "msg": run_.stderr.decode(), "out_fpath": Path(out_dump_fpath)}
    return results

 
def calculate_hetkmers(dump_fpath, out_fdir):
    name = dump_fpath.name
    out_fname = "{}_hetkmers".format(name)
    out_fpath = out_fdir/out_fname
    cmd = "smudgeplot.py hetkmers -o {} < {}"
    cmd = cmd.format(out_fpath, dump_fpath)
    sequence_file = dump_fpath.parent / "{}_sequences.tsv".format(str(out_fname))
    print(sequence_file)
    if sequence_file.exists():
        results = {"command": cmd, "returncode": 99, "name": name,
                   "msg": "output file exists already", 
                   "out_fpath": sequence_file}
    else:
        run_ = run(cmd, shell=True, capture_output=True)
        results = {"command": cmd, "returncode": run_.returncode,
                   "msg": run_.stderr.decode(), "out_fpath": str(sequence_file)}
    return results