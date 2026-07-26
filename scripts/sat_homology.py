#!/usr/bin/env python3
"""All-against-all homology between satDNA multimers, RepeatMasker-style.

Python 3 reimplementation of rm_homology.py (F. J. Ruiz-Ruano, satminer):
RepeatMasker is run once per sequence, using that single sequence as the
library and the whole multimer set as the query. Running it sequence by
sequence -- rather than one pass with the full library -- is the point: in a
single pass the library entries compete for each region and only the best hit
is reported, so the similarity between two related satDNAs is lost exactly
when it matters.

Output is a directional pairwise table (query = masked sequence, subject =
library sequence) with identity (100 - RepeatMasker divergence) and the
fraction of the shorter monomer covered by the alignment. sat_group.py turns
that into variants / families / superfamilies.

Usage:
  sat_homology.py --dimers dimers.fasta --lengths monomers.tsv \
      --out pairwise.tsv [--threads 8] [--workdir DIR] [--keep-workdir]
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor


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


def base_name(field):
    """'R1_CL1_1_1#Satellite' -> 'R1_CL1_1_1' (RepeatMasker keeps the class tag)."""
    return field.split("#", 1)[0]


def parse_rm_out(path):
    """Yield hits from a RepeatMasker .out file as dicts.

    Columns: SW  %div  %del  %ins  query  qbeg  qend  (qleft)  strand  repeat
             class/family  rbeg  rend  (rleft)  ID  [*]
    """
    if not os.path.isfile(path):
        return
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue                      # header / blank / continuation lines
            try:
                yield {
                    "sw": int(f[0]),
                    "div": float(f[1]),
                    "query": base_name(f[4]),
                    "qbeg": int(f[5]),
                    "qend": int(f[6]),
                    "subject": base_name(f[9]),
                }
            except ValueError:
                continue


def run_one(name, seq, dimers, workdir, engine, cutoff, extra):
    """RepeatMasker the whole multimer set with a single-sequence library."""
    jobdir = os.path.join(workdir, "rm_%s" % name.replace("/", "_"))
    os.makedirs(jobdir, exist_ok=True)
    lib = os.path.join(jobdir, "lib.fa")
    with open(lib, "w") as fh:                # '#Satellite' so RM reports a class
        fh.write(">%s#Satellite\n%s\n" % (name, seq))

    cmd = ["RepeatMasker", "-lib", lib, "-dir", jobdir, "-pa", "1",
           "-nolow", "-no_is", "-s", "-cutoff", str(cutoff)]
    if engine:
        cmd += ["-e", engine]
    cmd += extra + [dimers]

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True)
    out = os.path.join(jobdir, os.path.basename(dimers) + ".out")
    if proc.returncode != 0 and not os.path.isfile(out):
        sys.stderr.write("[sat_homology] RepeatMasker failed for %s:\n%s\n"
                         % (name, proc.stdout[-2000:]))
        return name, []
    return name, list(parse_rm_out(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dimers", required=True, help="multimer FASTA (from dimerator.py)")
    ap.add_argument("--lengths", required=True,
                    help="TSV with seq_id and monomer_len columns")
    ap.add_argument("--out", required=True, help="output pairwise TSV")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--engine", default=None, help="RepeatMasker -e (rmblast, hmmer, ...)")
    ap.add_argument("--cutoff", type=int, default=200,
                    help="RepeatMasker -cutoff Smith-Waterman score (default 200; "
                         "lower than RM's 225 default so short satDNAs are not lost)")
    ap.add_argument("--workdir", default=None, help="scratch dir (default: temp dir)")
    ap.add_argument("--keep-workdir", action="store_true")
    ap.add_argument("--rm-extra", default="",
                    help="extra RepeatMasker arguments, space separated")
    args = ap.parse_args()

    mono = {}
    with open(args.lengths) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                mono[row["seq_id"]] = int(row["monomer_len"])
            except (KeyError, ValueError, TypeError):
                continue

    seqs = [(base_name(h), s) for h, s in read_fasta(args.dimers)]
    fields = ["query", "subject", "identity", "divergence", "aln_len",
              "coverage", "sw_score"]

    if not seqs:
        with open(args.out, "w", newline="") as fh:
            csv.DictWriter(fh, fieldnames=fields, delimiter="\t").writeheader()
        sys.stderr.write("[sat_homology] empty input -> empty %s\n" % args.out)
        return

    workdir = args.workdir or tempfile.mkdtemp(prefix="sat_homology_")
    os.makedirs(workdir, exist_ok=True)
    dimers = os.path.abspath(args.dimers)
    extra = args.rm_extra.split() if args.rm_extra else []

    # Best directional hit per (query, subject) pair, ranked by identity then length.
    best = {}
    try:
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
            futures = [pool.submit(run_one, name, seq, dimers, workdir,
                                   args.engine, args.cutoff, extra)
                       for name, seq in seqs]
            for done, fut in enumerate(futures, 1):
                name, hits = fut.result()
                for h in hits:
                    q, s = h["query"], h["subject"]
                    if q == s:
                        continue              # self-hit, not informative
                    aln = abs(h["qend"] - h["qbeg"]) + 1
                    ref = min(mono.get(q, aln), mono.get(s, aln)) or aln
                    cov = min(1.0, float(aln) / ref)
                    rec = (round(100.0 - h["div"], 2), h["div"], aln,
                           round(cov, 3), h["sw"])
                    key = (q, s)
                    if key not in best or rec[0] > best[key][0] or (
                            rec[0] == best[key][0] and rec[2] > best[key][2]):
                        best[key] = rec
                if done % 25 == 0 or done == len(futures):
                    sys.stderr.write("[sat_homology] %d/%d sequences masked\n"
                                     % (done, len(futures)))
    finally:
        if not args.keep_workdir and args.workdir is None:
            shutil.rmtree(workdir, ignore_errors=True)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        for (q, s), (ident, div, aln, cov, sw) in sorted(best.items()):
            w.writerow([q, s, ident, div, aln, cov, sw])

    sys.stderr.write("[sat_homology] %d sequence(s), %d homologous pair(s) -> %s\n"
                     % (len(seqs), len(best), args.out))


if __name__ == "__main__":
    main()
