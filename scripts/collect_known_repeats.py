#!/usr/bin/env python3
"""Collect a representative "known repeats" FASTA from a RepeatExplorer2 run.

Python 3 port of satMiner's rexp_get_contigs_re2.py + extract_seq.py
(F. J. Ruiz-Ruano), with the BioPython dependency dropped and the cluster
directory passed explicitly instead of relying on the current working dir.

The 50%-cumulative-coverage selection (the actual satMiner logic) is kept
identical: for every cluster, contigs are ranked by coverage and the most
abundant ones are kept until their coverage sums to half of that cluster's
total. This yields a compact but representative set of the repeats RE2 has
already found, which we then use as a depletion reference.

Usage: collect_known_repeats.py <clusters_dir> <out.fasta>
       <clusters_dir> = .../re_output/seqclust/clustering/clusters
"""
import glob
import os
import re
import sys

clusters_dir, out_fasta = sys.argv[1], sys.argv[2]


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
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


# 1. Read every cluster's contigs.info.fasta into one {id: (header, seq)} table.
records = {}
info_files = sorted(glob.glob(os.path.join(clusters_dir, "dir_CL*", "contigs.info.fasta")))
for f in info_files:
    for header, seq in read_fasta(f):
        records[header.split()[0]] = (header, seq)

# 2. Parse cluster id, contig id and coverage from each header.
#    Headers look like "CL2Contig1 <n> <len> <coverage>". satMiner read the
#    coverage from a fixed positional field; we take the last numeric token,
#    which is where RE2 puts read depth and is robust to spacing differences.
cov_by_cluster = {}  # {CL: {contig: coverage}}
for cid, (header, _seq) in records.items():
    m = re.match(r"CL(\d+)Contig(\d+)", cid)
    if not m:
        continue
    cl, contig = m.group(1), m.group(2)
    cov = None
    for tok in reversed(header.split()):
        try:
            cov = float(tok)
            break
        except ValueError:
            continue
    if cov is None:
        continue
    cov_by_cluster.setdefault(cl, {})[contig] = cov

# 3. Per cluster keep top contigs until cumulative coverage >= 50% of the total.
selected = []
for cl, contigs in cov_by_cluster.items():
    ranked = sorted(contigs.items(), key=lambda kv: kv[1], reverse=True)
    half = sum(cov for _, cov in ranked) / 2.0
    acc = 0.0
    for contig, cov in ranked:
        if acc < half:                       # still below the halfway mark -> keep it
            selected.append("CL%sContig%s" % (cl, contig))
        acc += cov

# 4. Write the selected contig sequences out as the depletion reference.
n = 0
with open(out_fasta, "w") as w:
    for cid in selected:
        if cid in records:
            header, seq = records[cid]
            w.write(">%s\n%s\n" % (header, seq))
            n += 1

print("[collect_known_repeats] %d clusters, %d contigs selected -> %s"
      % (len(cov_by_cluster), n, out_fasta))
