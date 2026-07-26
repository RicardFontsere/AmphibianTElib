#!/usr/bin/env python3
"""Turn satDNA monomers into tandem multimers.

Python 3 reimplementation of dimerator.py (F. J. Ruiz-Ruano, satminer /
ngs-protocols) with one refinement taken from the more recent satellitome
papers: instead of always doubling the monomer, each monomer is repeated as
many times as needed to reach --min-len (at least twice).

Why this matters: RepeatMasker cannot anchor an alignment inside a monomer
shorter than its word/seed length, and a monomer is a circular permutation --
a read can start anywhere in it. A multimer both restores the junction
sequence and gives RepeatMasker enough length to align. A 21 bp monomer needs
ten copies to clear 200 bp; doubling it would leave it invisible.

Usage:
  dimerator.py --min-len 200 in.fasta out.fasta [--class Satellite] [--table t.tsv]
"""
import argparse
import csv
import math
import sys


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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infasta")
    ap.add_argument("outfasta")
    ap.add_argument("--min-len", type=int, default=200,
                    help="minimum multimer length in nt (default 200)")
    ap.add_argument("--min-copies", type=int, default=2,
                    help="minimum number of tandem copies (default 2)")
    ap.add_argument("--repeat-class", default=None,
                    help="append '#CLASS' to each name, e.g. Satellite, so that "
                         "RepeatMasker -lib reports a class in its .out file")
    ap.add_argument("--table", default=None, help="optional TSV log of copy numbers")
    ap.add_argument("--line-width", type=int, default=60)
    args = ap.parse_args()

    rows, n = [], 0
    with open(args.outfasta, "w") as out:
        for header, seq in read_fasta(args.infasta):
            name = header.split()[0]
            base = name.split("#", 1)[0]      # drop any class already present
            seq = seq.strip().upper()
            if not seq:
                continue
            copies = max(args.min_copies,
                         int(math.ceil(float(args.min_len) / len(seq))))
            multimer = seq * copies
            tag = "%s#%s" % (base, args.repeat_class) if args.repeat_class else name
            out.write(">%s\n" % tag)
            for i in range(0, len(multimer), args.line_width):
                out.write(multimer[i:i + args.line_width] + "\n")
            rows.append({"seq_id": base, "monomer_len": len(seq),
                         "copies": copies, "multimer_len": len(multimer)})
            n += 1

    if args.table:
        with open(args.table, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["seq_id", "monomer_len", "copies",
                                               "multimer_len"], delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    if rows:
        mx = max(r["copies"] for r in rows)
        sys.stderr.write("[dimerator] %d monomer(s) -> multimers >= %d nt "
                         "(max %d copies) -> %s\n" % (n, args.min_len, mx, args.outfasta))
    else:
        sys.stderr.write("[dimerator] no input sequences -> empty %s\n" % args.outfasta)


if __name__ == "__main__":
    main()
