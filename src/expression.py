
from csv import DictReader
from math import log10 as log


def calculate_sample_estimators(filepath, estimators, units, exclude, binary=False):
     with open(filepath) as fhand:
          raw_values = [float(row[units]) for row in DictReader(fhand, delimiter="\t") if row["Reference"] not in exclude]
          universe = len(raw_values)
          if binary:
               raw_values = [1 if float(value) >= 1 else 0 for value in raw_values]
          N = sum(raw_values)
          values = [(float(value)/N) * log(float(value)/N) if value > 0 else 0 for value in raw_values]
          diversity_value =  -sum(value for value in values if value != 0)
          specifity = log(universe) - diversity_value
          values = {"diversity": diversity_value, "specifity": specifity, "file": filepath}
          estimators[filepath.stem] = values
          return values, universe
