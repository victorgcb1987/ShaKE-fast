import math
import subprocess
import zlib

from dataclasses import dataclass

import numpy as np

from src.kmc import build_dump_argv
from src.kmer import LOG10_2


#Tunable knobs for hetkmer detection. Kept as plain constants rather than CLI
#flags, per the decision to keep the CLI surface to just --pooled.
POOLED_MAX_COMPONENT_SIZE = 6
NONPOOLED_MAX_COMPONENT_SIZE = 2
NONPOOLED_MIN_COVERAGE_RATIO = 0.2
#Caps the per-batch n x 3k candidate matrix in find_hamming1_edges. Unbatched
#at k=21, n=100M would need ~50GB - infeasible, so candidates are generated
#in chunks of this many k-mers at a time instead.
NEIGHBOR_BATCH_SIZE = 2_000_000

_BASE_TO_CODE = np.zeros(256, dtype=np.uint64)
for _i, _base in enumerate(b"ACGT"):
    _BASE_TO_CODE[_base] = _i
_CODE_TO_BASE = np.array(list("ACGT"))


@dataclass
class KmerTable:
    kmers: np.ndarray     # uint64, ascending sorted, unique, canonical-encoded (2 bits/base)
    counts: np.ndarray    # int64, aligned with kmers
    k: int
    kind: str = None
    merged: bool = False

    def universe_size(self):
        return int(self.kmers.shape[0])


#---------- 2-bit encoding ----------
#A k-mer is packed into a uint64 as k 2-bit base codes (A=0,C=1,G=2,T=3), most
#significant base first - this preserves lexicographic string ordering, which
#is what lets np.minimum(encoded, revcomp_encoded) reproduce KMC's own
#canonical-form rule (min(kmer, revcomp(kmer))) directly on the integers.

def encode_kmer_strings(kmer_bytes, k):
    n = len(kmer_bytes)
    if n == 0:
        return np.empty(0, dtype=np.uint64)
    flat = b"".join(kmer_bytes)
    arr = np.frombuffer(flat, dtype=np.uint8).reshape(n, k)
    codes2d = _BASE_TO_CODE[arr]
    encoded = np.zeros(n, dtype=np.uint64)
    for pos in range(k):
        encoded = (encoded << np.uint64(2)) | codes2d[:, pos]
    return encoded


def decode_kmer_codes(codes, k):
    codes = np.asarray(codes, dtype=np.uint64)
    out = np.empty((codes.shape[0], k), dtype="<U1")
    remaining = codes.copy()
    for pos in range(k - 1, -1, -1):
        idx = (remaining & np.uint64(3)).astype(np.int64)
        out[:, pos] = _CODE_TO_BASE[idx]
        remaining = remaining >> np.uint64(2)
    return ["".join(row) for row in out]


def revcomp_encoded(codes, k):
    full_mask = np.uint64((1 << (2 * k)) - 1)
    comp = codes ^ full_mask
    out = np.zeros_like(comp)
    for pos in range(k):
        group = (comp >> np.uint64(2 * pos)) & np.uint64(3)
        out |= group << np.uint64(2 * (k - 1 - pos))
    return out


def canonicalize_encoded(codes, k):
    return np.minimum(codes, revcomp_encoded(codes, k))


#---------- loading (replaces writing+reading a .dump file entirely) ----------

def load_kmer_table(db_fpath, k, threads=6, lower_bound=1, upper_bound=9999999999, merged=False):
    argv = build_dump_argv(db_fpath, threads, lower_bound, upper_bound)
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    raw_stdout, raw_stderr = proc.communicate()
    results = {"command": " ".join(argv), "returncode": proc.returncode,
               "msg": raw_stderr.decode(), "name": str(db_fpath)}
    if proc.returncode != 0:
        empty = KmerTable(np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.int64), k, merged=merged)
        return empty, results
    table = parse_dump_bytes(raw_stdout, k, merged=merged)
    #kmc_tools is only trusted for -ci (see build_dump_argv): enforce both
    #bounds here so the result never depends on how kmc_tools parses them.
    keep = (table.counts >= lower_bound) & (table.counts <= upper_bound)
    if not keep.all():
        table = KmerTable(table.kmers[keep], table.counts[keep], k, merged=merged)
    return table, results


def parse_dump_bytes(raw, k, merged=False):
    if not raw:
        return KmerTable(np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.int64), k, merged=merged)
    lines = raw.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    kmer_bytes = [None] * len(lines)
    counts = np.empty(len(lines), dtype=np.int64)
    for i, line in enumerate(lines):
        kmer, count = line.split(b"\t")
        kmer_bytes[i] = kmer
        counts[i] = int(count)
    codes = encode_kmer_strings(kmer_bytes, k)
    if codes.shape[0] > 1 and not np.all(codes[:-1] <= codes[1:]):
        #Defensive: `dump -s` should already yield ascending order, but don't
        #trust it blindly.
        order = np.argsort(codes, kind="stable")
        codes, counts = codes[order], counts[order]
    return KmerTable(codes, counts, k, merged=merged)


#---------- universe size ----------

def union_universe_size(tables):
    tables = [table for table in tables if table.kmers.shape[0] > 0]
    if not tables:
        return 0
    if len(tables) == 1:
        return tables[0].universe_size()
    combined = np.concatenate([table.kmers for table in tables])
    return int(np.unique(combined).shape[0])


#---------- Hamming-distance-1 hetkmer detection ----------

def hamming1_neighbor_masks(k):
    pos_shift = (2 * (k - 1 - np.arange(k))).astype(np.uint64)
    deltas = np.array([1, 2, 3], dtype=np.uint64)
    return (deltas[None, :] << pos_shift[:, None]).reshape(-1)


def find_hamming1_edges(table, batch_size=NEIGHBOR_BATCH_SIZE):
    n = table.kmers.shape[0]
    if n == 0:
        return np.empty((0, 2), dtype=np.int64)
    masks = hamming1_neighbor_masks(table.k)
    n_masks = masks.shape[0]
    batches = []
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        batch = table.kmers[start:end]
        candidates = batch[:, None] ^ masks[None, :]
        canon = canonicalize_encoded(candidates.reshape(-1), table.k)
        idx = np.clip(np.searchsorted(table.kmers, canon), 0, n - 1)
        found = table.kmers[idx] == canon
        b = end - start
        rows_local, cols = np.nonzero(found.reshape(b, n_masks))
        matched = idx.reshape(b, n_masks)[rows_local, cols]
        rows_global = rows_local + start
        a_idx = np.minimum(rows_global, matched)
        b_idx = np.maximum(rows_global, matched)
        valid = a_idx != b_idx
        if valid.any():
            batches.append(np.stack([a_idx[valid], b_idx[valid]], axis=1))
    if not batches:
        return np.empty((0, 2), dtype=np.int64)
    edges = np.concatenate(batches, axis=0)
    return np.unique(edges, axis=0)


def _find(parent, x):
    root = x
    while parent[root] != root:
        root = parent[root]
    while parent[x] != root:
        parent[x], x = root, parent[x]
    return root


def build_components(table_size, edges):
    #Purpose-built array union-find local to this module: only the k-mers
    #actually touched by a Hamming-1 edge ever need to enter it, so this
    #stays cheap regardless of the total k-mer count. Kept separate from
    #src/utils.py:UnionFind (dict-based/recursive, still used as-is by
    #legacy_scripts/group_kmers_by_hetkmers.py).
    if edges.shape[0] == 0:
        return {}
    involved = np.unique(edges)
    local = np.searchsorted(involved, edges)
    m = involved.shape[0]
    parent = np.arange(m, dtype=np.int64)
    for a, b in local:
        ra, rb = _find(parent, a), _find(parent, b)
        if ra != rb:
            parent[rb] = ra
    for i in range(m):
        parent[i] = _find(parent, i)
    groups = {}
    for local_idx, root in enumerate(parent):
        groups.setdefault(int(root), []).append(int(involved[local_idx]))
    return groups


def merge_hetkmers_native(table, pooled):
    max_size = POOLED_MAX_COMPONENT_SIZE if pooled else NONPOOLED_MAX_COMPONENT_SIZE
    edges = find_hamming1_edges(table)
    groups = build_components(table.kmers.shape[0], edges)
    to_drop = []
    add_kmers = []
    add_counts = []
    kept = 0
    rejected = 0
    for members in groups.values():
        if len(members) > max_size:
            rejected += 1
            continue
        if not pooled:
            #max_size == 2 here, so members has exactly 2 entries.
            a, b = members
            count_a, count_b = int(table.counts[a]), int(table.counts[b])
            low, high = min(count_a, count_b), max(count_a, count_b)
            if high == 0 or (low / high) < NONPOOLED_MIN_COVERAGE_RATIO:
                rejected += 1
                continue
        to_drop.extend(members)
        add_kmers.append(int(table.kmers[members].min()))
        add_counts.append(int(table.counts[members].sum()))
        kept += 1
    stats = {"components_found": len(groups), "components_merged": kept, "components_rejected": rejected}
    if not to_drop:
        return table, stats
    keep_mask = np.ones(table.kmers.shape[0], dtype=bool)
    keep_mask[to_drop] = False
    new_kmers = np.concatenate([table.kmers[keep_mask], np.array(add_kmers, dtype=np.uint64)])
    new_counts = np.concatenate([table.counts[keep_mask], np.array(add_counts, dtype=np.int64)])
    order = np.argsort(new_kmers, kind="stable")
    merged_table = KmerTable(new_kmers[order], new_counts[order], table.k,
                              kind=table.kind, merged=table.merged)
    return merged_table, stats


#---------- Shannon diversity ----------

def shannon_diversity_log10(counts, presence=False):
    #The part of the Shannon estimators that only depends on the sample's own
    #counts - not on the universe size.
    values = (counts >= 1).astype(np.float64) if presence else counts.astype(np.float64)
    N = values.sum()
    p = values[values > 0] / N
    return -float(np.sum(p * np.log10(p)))


def shannon_from_diversity_log10(diversity_log10, universe_size):
    #Completes the estimators once the universe size is known.
    specifity_log10 = math.log10(universe_size) - diversity_log10
    diversity_log2 = diversity_log10 / LOG10_2
    specifity_log2 = math.log2(universe_size) - diversity_log2
    return {"diversity_log10": diversity_log10, "specifity_log10": specifity_log10,
            "diversity_log2": diversity_log2, "specifity_log2": specifity_log2}


def evenness_from_diversity_log10(diversity_log10, n_positive, universe_size):
    #Shannon diversity rescaled to 0-1 (the ratio is the same in any log base):
    #  pielou_evenness   = H / log(n)  - n: kmers/genes actually present in the
    #                      sample, so it only depends on the sample itself;
    #  universe_evenness = H / log(U)  - U: the universe size, so it also
    #                      reflects how much of the universe the sample covers.
    #NaN when the denominator is 0 (n <= 1 or U <= 1). Not clipped: H/log(U)
    #can marginally exceed 1 if a sample has more kmers than its universe
    #(possible after het-kmer merging).
    nan = float("nan")
    return {"pielou_evenness": diversity_log10 / math.log10(n_positive) if n_positive > 1 else nan,
            "universe_evenness": diversity_log10 / math.log10(universe_size) if universe_size > 1 else nan}


def shannon_estimators_from_counts(counts, universe_size, presence=False):
    return shannon_from_diversity_log10(shannon_diversity_log10(counts, presence), universe_size)


#---------- Kolmogorov complexity ratio, fully in-memory ----------

def counts_to_binary_bytes(values, bit_width, num_zeros=0):
    n = values.shape[0]
    out = np.zeros((n, bit_width + 1), dtype=np.uint8)
    out[:, -1] = 10  # '\n'
    v = values.astype(np.uint64)
    for bit in range(bit_width):
        shift = bit_width - 1 - bit
        out[:, bit] = 48 + ((v >> np.uint64(shift)) & np.uint64(1)).astype(np.uint8)
    payload = out.tobytes()
    if num_zeros > 0:
        payload += (b"0" * bit_width + b"\n") * num_zeros
    return payload


def kolmogorov_ratio(raw_bytes, compresslevel=6):
    if not raw_bytes:
        return 0.0
    return len(zlib.compress(raw_bytes, compresslevel)) / len(raw_bytes)


def kolmogorov_from_table(table, universe_size, presence=False):
    n = table.counts.shape[0]
    num_zeros = max(universe_size - n, 0)
    if presence:
        payload = b"01\n" * n + b"00\n" * num_zeros
    else:
        payload = counts_to_binary_bytes(table.counts, 30, num_zeros)
    return kolmogorov_ratio(payload)


class StreamingKolmogorov:
    """kolmogorov_from_table split in two so the table doesn't have to be kept.

    The sample's own counts are compressed as soon as it is built (in chunks,
    never materializing the whole payload); only the compressor state (a few
    hundred KB) is kept. When the universe size is finally known, `finish`
    feeds the missing zeros (also in chunks) and returns the same ratio
    kolmogorov_from_table gives.
    """
    BIT_WIDTH = 30
    CHUNK = 262144

    def __init__(self, counts=None, presence=False, compresslevel=6):
        self.presence = presence
        self.n = 0
        self._compressor = zlib.compressobj(compresslevel)
        self._compressed = 0
        self._raw = 0
        if counts is not None:
            for start in range(0, counts.shape[0], self.CHUNK):
                self.feed_counts(counts[start:start + self.CHUNK])

    def feed_counts(self, part):
        self.n += int(part.shape[0])
        if self.presence:
            self._feed(b"01\n" * part.shape[0])
        else:
            self._feed(counts_to_binary_bytes(part, self.BIT_WIDTH))

    def _feed(self, payload):
        self._raw += len(payload)
        self._compressed += len(self._compressor.compress(payload))

    def finish_sizes(self, universe_size):
        #Pads with zeros up to the universe size and returns
        #(compressed_size, raw_size). Can only be called once.
        zero_line = b"00\n" if self.presence else b"0" * self.BIT_WIDTH + b"\n"
        remaining = max(universe_size - self.n, 0)
        while remaining > 0:
            step = min(remaining, self.CHUNK)
            self._feed(zero_line * step)
            remaining -= step
        self._compressed += len(self._compressor.flush())
        return self._compressed, self._raw

    def finish(self, universe_size):
        compressed, raw = self.finish_sizes(universe_size)
        return compressed / raw if raw else 0.0


#---------- normalized Kolmogorov ----------
#The raw ratio lives in a very narrow band (an even profile can't go below
#deflate's ~1000:1 cap, a random 30-bit one can't go above ~0.19 of the count
#lines) and depends on depth and universe size. The normalized version
#places the sample between two references with the same n and universe size:
#  floor   - the presence/absence payload (every count the same);
#  ceiling - i.i.d. geometric counts with the sample's own mean, the maximum
#            entropy distribution over non-negative integers at a given mean.
REFERENCE_SEED = 20240601


def reference_compressed_size(n, mean_count, universe_size, seed=REFERENCE_SEED):
    #Compressed size of the geometric reference. Deterministic (fixed seed,
    #fixed chunking), so the same inputs always give the same number.
    reference = StreamingKolmogorov(presence=False)
    if n > 0:
        rng = np.random.default_rng(seed)
        p = min(1.0, 1.0 / max(mean_count, 1.0))
        cap = (1 << StreamingKolmogorov.BIT_WIDTH) - 1
        done = 0
        while done < n:
            step = min(n - done, StreamingKolmogorov.CHUNK)
            reference.feed_counts(np.minimum(rng.geometric(p, size=step), cap).astype(np.int64))
            done += step
    return reference.finish_sizes(universe_size)[0]


def normalize_kolmogorov(sample_size, floor_size, reference_size):
    #(C - floor) / (reference - floor). NaN when it is undefined (no k-mers,
    #or a reference that isn't more complex than the floor). Not clipped: a
    #sample can come out slightly above 1 because deflate doesn't reach the
    #true entropy of the reference.
    denominator = reference_size - floor_size
    if denominator <= 0:
        return float("nan")
    return (sample_size - floor_size) / denominator


def kolmogorov_estimators_from_counts(counts, universe_size):
    #kolmogorov, kolmogorov_presence and kolmogorov_norm in one pass over the
    #counts (same numbers as kolmogorov_from_table for the first two).
    regular = StreamingKolmogorov(counts, presence=False)
    presence = StreamingKolmogorov(counts, presence=True)
    return kolmogorov_estimators_from_streams(regular, presence, counts.sum(), universe_size)


def kolmogorov_estimators_from_streams(regular, presence, total_count, universe_size):
    n = regular.n
    sample_size, raw = regular.finish_sizes(universe_size)
    floor_size, raw_presence = presence.finish_sizes(universe_size)
    if n == 0:
        norm = float("nan")
    else:
        reference_size = reference_compressed_size(n, total_count / n, universe_size)
        norm = normalize_kolmogorov(sample_size, floor_size, reference_size)
    return {"kolmogorov": sample_size / raw if raw else 0.0,
            "kolmogorov_presence": floor_size / raw_presence if raw_presence else 0.0,
            "kolmogorov_norm": norm}
