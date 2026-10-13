
from csv import DictReader
from math import log10, log2

import numpy as np

from src.kmerdata import CountHistogram, count_complexity_estimators


def calculate_sample_estimators(filepath, estimators, units, exclude, binary=False):
     with open(filepath) as fhand:
          raw_values = [float(row[units]) for row in DictReader(fhand, delimiter="\t") if row["Reference"] not in exclude]
          universe = len(raw_values)
          if binary:
               raw_values = [1 if float(value) >= 1 else 0 for value in raw_values]
          n_positive = sum(1 for value in raw_values if value > 0)
          N = sum(raw_values)
          values_log10 = [(float(value)/N) * log10(float(value)/N) if value > 0 else 0 for value in raw_values]
          diversity_log10 =  -sum(value for value in values_log10 if value != 0)
          specifity_log10 = log10(universe) - diversity_log10
          values_log2 = [(float(value)/N) * log2(float(value)/N) if value > 0 else 0 for value in raw_values]
          diversity_log2 =  -sum(value for value in values_log2 if value != 0)
          specifity_log2 = log2(universe) - diversity_log2
          values = {"diversity": diversity_log10, "specifity": specifity_log10,
                    "diversity_log10": diversity_log10, "specifity_log10": specifity_log10,
                    "diversity_log2": diversity_log2, "specifity_log2": specifity_log2,
                    "n_positive": n_positive,
                    "file": filepath}
          estimators[filepath.stem] = values
          return values, universe


def count_complexity_from_expression_file(filepath, units, exclude):
     #Count complexity of an expression table: the values (rounded to 3 decimals
     #and scaled to integers, as the old Kolmogorov encoding did) play the role
     #of the counts and the genes of the table are the universe (zero values
     #are the absent elements).
     with open(filepath) as fhand:
          values = [round(float(row[units]), 3) * 1000 for row in DictReader(fhand, delimiter="\t")
                    if row["Reference"] not in exclude]
     integers = np.trunc(np.array(values, dtype=np.float64)).astype(np.int64)
     return count_complexity_estimators(CountHistogram(integers), len(values))
