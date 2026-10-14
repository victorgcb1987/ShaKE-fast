#!/usr/bin/env python
"""Unit checks of the count complexity functions in src/kmerdata.py (no kmc needed):

  * the decomposition L/U = H2(n/U) + (n/U) * H_c equals the direct sum over the histogram;
  * closed-form H_max(mean) equals the entropy of the geometric distribution (numerically);
  * spectrum_entropy_norm is about 1 for geometric counts and 0 for constant counts (NaN: no spread);
  * zeros: padding the same sample to a larger universe lowers the presence part as expected;
  * edge cases: empty sample, one count value, universe smaller than n, chunked histogram.

Usage:
    python validation/check_spectrum_entropy.py
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.kmerdata import (CountHistogram, spectrum_entropy_estimators, max_entropy_nonnegative)

failures = []


def check(name, condition, detail=""):
    print("{} {}{}".format("ok  " if condition else "FAIL", name, " ({})".format(detail) if detail else ""))
    if not condition:
        failures.append(name)


rng = np.random.default_rng(3)
# 1. decomposition against the direct sum
counts = rng.geometric(0.15, size=50_000) + 1          # minimum count 2, like a lower bound of 2
for universe in (50_000, 80_000, 500_000):
    est = spectrum_entropy_estimators(CountHistogram(counts), universe)
    values, freq = np.unique(counts, return_counts=True)
    direct = float(np.sum(freq * np.log2(universe / freq)))
    zeros = universe - counts.shape[0]
    if zeros:
        direct += zeros * math.log2(universe / zeros)
    check("decomposition, U = {}".format(universe), abs(est["spectrum_entropy_universe"] - direct / universe) < 1e-12,
          "{:.6f}".format(est["spectrum_entropy_universe"]))
# 2. closed form of the maximum entropy
for mu in (0.5, 3.0, 25.0):
    k = np.arange(0, 20000)
    p = (1 / (1 + mu)) * (mu / (1 + mu)) ** k
    numeric = float(-np.sum(p[p > 0] * np.log2(p[p > 0])))
    check("H_max closed form, mean {}".format(mu), abs(max_entropy_nonnegative(mu) - numeric) < 1e-9)
# 3. normalised value
geometric = rng.geometric(1 / 8, size=2_000_000)
norm = spectrum_entropy_estimators(CountHistogram(geometric), geometric.shape[0])["spectrum_entropy_norm"]
check("norm about 1 for geometric counts", abs(norm - 1) < 0.01, "{:.4f}".format(norm))
constant = np.full(1000, 7)
check("norm is NaN for constant counts", math.isnan(spectrum_entropy_estimators(CountHistogram(constant), 1000)["spectrum_entropy_norm"]))
# 4. zeros and universe
est_small = spectrum_entropy_estimators(CountHistogram(counts), counts.shape[0])
est_big = spectrum_entropy_estimators(CountHistogram(counts), counts.shape[0] * 4)
check("own universe: presence part is 0", est_small["occupancy_entropy"] == 0.0)
check("shared universe: presence part = H2(1/4)", abs(est_big["occupancy_entropy"] - (-(0.25 * math.log2(0.25)) - 0.75 * math.log2(0.75))) < 1e-12)
check("count entropy and norm do not depend on the universe",
      est_small["spectrum_entropy"] == est_big["spectrum_entropy"] and est_small["spectrum_entropy_norm"] == est_big["spectrum_entropy_norm"])
# 5. edge cases
empty = spectrum_entropy_estimators(CountHistogram(np.array([], dtype=np.int64)), 100)
check("empty sample", empty["spectrum_entropy_universe"] == 0.0 and math.isnan(empty["spectrum_entropy_norm"]))
clamped = spectrum_entropy_estimators(CountHistogram(counts), 10)
check("universe smaller than n is clamped to n", clamped["occupancy_entropy"] == 0.0)
zeros_in_input = spectrum_entropy_estimators(CountHistogram(np.concatenate([counts, np.zeros(1000, dtype=counts.dtype)])), counts.shape[0] + 1000)
check("zeros in the counts are treated as absent", abs(zeros_in_input["spectrum_entropy"] - est_small["spectrum_entropy"]) < 1e-12)
chunked = CountHistogram()
chunked.CHUNK = 1000
chunked.feed_counts(counts)
check("chunked histogram equals one-shot histogram", chunked.hist == CountHistogram(counts).hist)

print("\n{} failures".format(len(failures)))
sys.exit(1 if failures else 0)
