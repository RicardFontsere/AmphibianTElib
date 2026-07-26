#!/usr/bin/env python3
"""Collect putative satDNA monomers from one or more RepeatExplorer2 runs.

TAREAN writes its consensus monomers to
    <re_output>/TAREAN_consensus_rank_<N>.fasta
with
    rank 1 = high confidence satellite
    rank 2 = low confidence satellite
    rank 3 = putative LTR element
    rank 4 = other (rDNA, low complexity, ...)
Only the ranks given with --ranks are kept (1,2 by default; ranks 3-4 are mobile
elements or rDNA, not satDNA). If those files are absent (older or rearranged
output trees) the per-cluster tarean/ directories are used as a fallback and the
confidence is reported as "unknown".

Cluster sizes are read best-effort from CLUSTER_TABLE.csv; they are informative
only -- the abundance that ends up in the catalogue is re-estimated later by
RepeatMasking a random read sample against the library.

Usage:
  collect_tarean_satellites.py --out sat.fasta --table sat.tsv \
      --ranks 1,2 R1=<re_output> R2=<re_output_r2>
"""
import argparse
import csv
import glob
import os
import re
import sys

RANK_LABEL = {
    "1": "high_confidence_satellite",
    "2": "low_confidence_satellite",
    "3": "putative_LTR",
    "4": "other",
}


def read_fasta(path):
    """Yield (header, sequence) pairs from a FASTA file (header without '>')."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def parse_cluster_table(re_dir):
    """Best-effort {cluster_id: {size, proportion, annotation}} from CLUSTER_TABLE.csv.

    The file carries a free-text preamble before the real header row, and RE2 has
    shipped both tab- and comma-separated variants, so the header row is located
    by content and the delimiter sniffed from it.
    """
    path = os.path.join(re_dir, "CLUSTER_TABLE.csv")
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        lines = fh.read().splitlines()

    head_i, delim = None, None
    for i, line in enumerate(lines):
        for d in ("\t", ","):
            cells = [c.strip().strip('"') for c in line.split(d)]
            low = [c.lower() for c in cells]
            if "cluster" in low and any(c.startswith("size") for c in low):
                head_i, delim = i, d
                break
        if head_i is not None:
            break
    if head_i is None:
        return {}

    header = [c.strip().strip('"') for c in lines[head_i].split(delim)]
    low = [c.lower() for c in header]

    def col(*names):
        for n in names:
            for j, c in enumerate(low):
                if c == n or c.startswith(n):
                    return j
        return None

    i_cl = col("cluster")
    i_size = col("size_adjusted", "size")
    i_prop = col("genome_proportion", "proportion")
    i_ann = col("tarean_annotation", "automatic_annotation", "annotation")

    out = {}
    for line in lines[head_i + 1:]:
        cells = [c.strip().strip('"') for c in line.split(delim)]
        if len(cells) <= i_cl or not cells[i_cl]:
            continue
        m = re.search(r"(\d+)", cells[i_cl])
        if not m:
            continue
        rec = {"size": "NA", "prop": "NA", "annotation": "NA"}
        for key, idx in (("size", i_size), ("prop", i_prop), ("annotation", i_ann)):
            if idx is not None and len(cells) > idx and cells[idx]:
                rec[key] = cells[idx]
        out[str(int(m.group(1)))] = rec
    return out


def collect_round(label, re_dir, ranks, min_len):
    """Return a list of monomer records for one RepeatExplorer2 output directory."""
    table = parse_cluster_table(re_dir)
    records, seen_per_cluster = [], {}

    rank_files = []
    for rank in ranks:
        p = os.path.join(re_dir, "TAREAN_consensus_rank_%s.fasta" % rank)
        if os.path.isfile(p):
            rank_files.append((rank, p))

    if not rank_files:
        # Fallback: per-cluster TAREAN consensus files, confidence unknown.
        fallback = sorted(glob.glob(os.path.join(
            re_dir, "seqclust", "clustering", "clusters", "dir_CL*", "tarean",
            "*consensus*.fasta")))
        if not fallback:
            sys.stderr.write(
                "[collect_tarean] WARNING: no TAREAN output found under %s\n" % re_dir)
            return records
        sys.stderr.write(
            "[collect_tarean] WARNING: no TAREAN_consensus_rank_*.fasta in %s, "
            "falling back to %d per-cluster consensus file(s); ranks not filtered\n"
            % (re_dir, len(fallback)))
        rank_files = [("NA", p) for p in fallback]

    for rank, path in rank_files:
        for header, seq in read_fasta(path):
            seq = re.sub(r"[^ACGTNacgtn]", "", seq).upper()
            if len(seq) < min_len:
                continue
            m = re.search(r"CL(\d+)", header)
            cluster = str(int(m.group(1))) if m else "NA"
            n = seen_per_cluster.get((cluster, rank), 0) + 1
            seen_per_cluster[(cluster, rank)] = n
            info = table.get(cluster, {})
            records.append({
                "seq_id": "%s_CL%s_%s_%d" % (label, cluster, rank, n),
                "round": label,
                "cluster": cluster,
                "tarean_rank": rank,
                "confidence": RANK_LABEL.get(rank, "unknown"),
                "monomer_len": len(seq),
                "re_cluster_size": info.get("size", "NA"),
                "re_cluster_prop_pct": info.get("prop", "NA"),
                "re_annotation": info.get("annotation", "NA"),
                "source": os.path.relpath(path, os.path.dirname(re_dir.rstrip("/"))),
                "original_header": header,
                "seq": seq,
            })
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rounds", nargs="+", metavar="LABEL=DIR",
                    help="RepeatExplorer2 output dirs, e.g. R1=re_output R2=re_output_r2")
    ap.add_argument("--out", required=True, help="output monomer FASTA")
    ap.add_argument("--table", required=True, help="output metadata TSV")
    ap.add_argument("--ranks", default="1,2",
                    help="comma-separated TAREAN ranks to keep (default 1,2)")
    ap.add_argument("--min-len", type=int, default=10,
                    help="skip monomers shorter than this (default 10)")
    args = ap.parse_args()

    ranks = [r.strip() for r in args.ranks.split(",") if r.strip()]

    records = []
    for spec in args.rounds:
        if "=" not in spec:
            ap.error("round argument must look like LABEL=DIR, got %r" % spec)
        label, re_dir = spec.split("=", 1)
        if not os.path.isdir(re_dir):
            sys.stderr.write("[collect_tarean] WARNING: missing dir %s\n" % re_dir)
            continue
        got = collect_round(label, re_dir, ranks, args.min_len)
        sys.stderr.write("[collect_tarean] %s: %d monomer(s) from %s\n"
                         % (label, len(got), re_dir))
        records.extend(got)

    fields = ["seq_id", "round", "cluster", "tarean_rank", "confidence", "monomer_len",
              "re_cluster_size", "re_cluster_prop_pct", "re_annotation", "source",
              "original_header"]
    with open(args.out, "w") as fa, open(args.table, "w", newline="") as tsv:
        w = csv.DictWriter(tsv, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in records:
            fa.write(">%s\n%s\n" % (r["seq_id"], r["seq"]))
            w.writerow(r)

    sys.stderr.write("[collect_tarean] %d putative satDNA monomer(s) -> %s\n"
                     % (len(records), args.out))
    if not records:
        sys.stderr.write("[collect_tarean] NOTE: empty library; downstream rules "
                         "will produce empty outputs for this species\n")


if __name__ == "__main__":
    main()
