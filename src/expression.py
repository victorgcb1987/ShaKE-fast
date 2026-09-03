
from csv import DictReader
from math import log10, log2


def calculate_sample_estimators(filepath, estimators, units, exclude, binary=False):
     with open(filepath) as fhand:
          raw_values = [float(row[units]) for row in DictReader(fhand, delimiter="\t") if row["Reference"] not in exclude]
          universe = len(raw_values)
          if binary:
               raw_values = [1 if float(value) >= 1 else 0 for value in raw_values]
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
                    "file": filepath}
          estimators[filepath.stem] = values
          return values, universe
