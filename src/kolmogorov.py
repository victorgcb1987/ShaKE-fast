import gzip
import os

from csv import DictReader
from pathlib import Path
from src.utils import get_universe_size
from subprocess import run


'''BINARY_LENGTH=30(a thousand of millions)!!!!
    Steps (expression):
        1.-Get Normalized data (TPM)
        2.-Round data and convert each value to binary
        3.-Compress data
        4.-Compressed size/original size
        
    Steps (kmer_counts)
        1.-Get universe size
        2.-Calculate difference between universe size and number of kmers in your sample
        3.-Get kmer counts and add 0s equal to this difference
        4.-Convert values to binary (value/total)
        5.-Compress data
        6.-Compressed size/original size'''


def get_universe_size_difference(filepath, universe_size):
    return universe_size - get_universe_size([filepath])


def convert_to_binary(number, presence=False):
    if presence:
        return format(1, '02b')
    return format(int(number), '030b')


# def create_kmer_binary_file(in_filepath, out_filepath, num_zeros):
#     compressed = "{}.gz".format(out_filepath)
#     with gzip.open(compressed, 'wb') as compressed_fhand:
#         with open(out_filepath, "w") as not_compressed_fhand:
#             with open(in_filepath) as in_fhand:
#                 generator = (convert_to_binary(line.rstrip().split()[1]) for line in in_fhand)
#                 for gen in generator:
#                     compressed_fhand.write(gen.encode()+b"\n")
#                     compressed_fhand.flush()
#                     not_compressed_fhand.write(gen+"\n")
#                     not_compressed_fhand.flush()
#                 for zero in range(num_zeros):
#                     compressed_fhand.write(format(0, '030b').encode()+b"\n")
#                     compressed_fhand.flush()
#                     not_compressed_fhand.write(format(0, '030b')+"\n")
#                     not_compressed_fhand.flush()
#     return compressed


def create_kmer_binary_file(in_filepath, out_filepath, num_zeros, presence=False):
    uncompressed_exists = Path(out_filepath).exists()
    if not uncompressed_exists:
        # NOTE: this used to call fhand.flush() after writing every single
        # k-mer line (and every padding zero). With tens/hundreds of
        # millions of k-mers that turns buffered sequential I/O into one
        # syscall per line and dominated the whole pipeline's runtime.
        # Python's buffered file object already flushes on close (end of
        # the `with` block), so no manual flushing is needed here.
        zero_line = (format(0, '02b') if presence else format(0, '030b')) + "\n"
        with open(out_filepath, "w") as not_compressed_fhand:
            with open(in_filepath) as in_fhand:
                lines = (convert_to_binary(line.rstrip().split()[1], presence=presence) + "\n" for line in in_fhand)
                not_compressed_fhand.writelines(lines)
                if num_zeros:
                    not_compressed_fhand.writelines(zero_line for _ in range(num_zeros))
    compressed = "{}.gz".format(out_filepath)
    compressed_exists = Path(compressed).exists()
    cmd = "gzip -c {} > {}".format(out_filepath, compressed)
    if not compressed_exists:
        run(cmd, shell=True)
    already_done = uncompressed_exists and compressed_exists
    return {"command": cmd, "returncode": 99 if already_done else 0,
            "msg": "output file exists already" if already_done else "",
            "out_fpath": compressed}


def calculate_kolmogorov(filepath_a, filepath_b):
    return float(os.stat(filepath_a).st_size/ os.stat(filepath_b).st_size)


def create_expression_binary_file(in_filepath, units, exclude, out_fpath, presence=False):
    # Same fix as create_kmer_binary_file: no per-line flush() (here it was
    # even worse - two flushes per row, one per output file), and write in
    # one batched call instead of one write() syscall per row.
    compressed = "{}.gz".format(out_fpath)
    cmd = "gzip-encode expression binary {} -> {}".format(out_fpath, compressed)
    if Path(compressed).exists():
        return {"command": cmd, "returncode": 99,
                "msg": "output file exists already", "out_fpath": compressed}
    with gzip.open(compressed, 'wb') as compressed_fhand:
        with open(out_fpath, "w") as not_compressed_fhand:
            with open(in_filepath) as fhand:
                lines = [convert_to_binary(round(float(row[units]), 3)*1000, presence=presence)
                         for row in DictReader(fhand, delimiter="\t") if row["Reference"] not in exclude]
                text = "\n".join(lines) + "\n" if lines else ""
                not_compressed_fhand.write(text)
                compressed_fhand.write(text.encode())
    return {"command": cmd, "returncode": 0, "msg": "", "out_fpath": compressed}


def calculate_kolmogorov_estimator(filepath, universe_size, estimators, group=None,
                                   sub=None, name=None, kind=None, units="TPM", presence=False,
                                   key="kolmogorov"):
    if presence:
        binary = "{}.presence.binary".format(str(filepath))
    else:
        binary = "{}.binary".format(str(filepath))
    if kind != "expression":
        num_zeros = get_universe_size_difference(filepath, universe_size)
        binary_results = create_kmer_binary_file(filepath, binary, num_zeros, presence=presence)
    else:
        binary_results = create_expression_binary_file(filepath, units, [], binary, presence=presence)
    compressed_file = binary_results["out_fpath"]
    kolmo = calculate_kolmogorov(compressed_file, binary)
    estimators[group][sub][name][key] = kolmo
    return {"command": binary_results["command"], "returncode": binary_results["returncode"],
            "msg": binary_results["msg"], "name": name, "out_fpath": compressed_file}


    