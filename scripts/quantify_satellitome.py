#!/usr/bin/env python3
"""Re-estimate satDNA abundance by RepeatMasking random read samples.

Cluster size in RepeatExplorer is not a comparable abundance estimate: it
depends on the read set, on the clustering threshold and on how the reads of
one family split over several clusters, and round-2 clusters are computed on
depleted reads. The published satellitomes therefore re-estimate abundance by
RepeatMasking a fixed random sample of reads against the multimerised library
and dividing masked bp by sampled bp. That is what this script summarises.

Reading is streaming and assumes RepeatMasker .out files are grouped by query
(they are): intervals are merged per (read, library entry) so overlapping hits
of one entry are not counted twice, and merged across entries for the total
satDNA fraction.

Usage:
  quantify_satellitome.py --totals totals.tsv --groups groups.tsv \
      --out abundance.tsv --summary summary.tsv rep1=rep1.fa.out rep2=rep2.fa.out
"""
import argparse
import csv
import math
import os
import sys


def base_name(field):
    """'RteSat01-176#Satellite' -> 'RteSat01-176'."""
    return field.split("#", 1)[0]


def merge(intervals):
    """Total length covered by a list of (start, end) inclusive intervals."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None:
            cur_s, cur_e = s, e
        elif s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def scan(path):
    """Return (masked_bp, hits, div_weight) per entry plus total masked bp."""
    per_entry = {}          # entry -> [masked_bp, hits, sum(div*bp)]
    total_masked = 0
    cur_query, buf = None, {}

    def flush():
        nonlocal total_masked
        if not buf:
            return
        allint = []
        for entry, ivs in buf.items():
            bp = merge([(s, e) for s, e, _d in ivs])
            rec = per_entry.setdefault(entry, [0, 0, 0.0])
            rec[0] += bp
            rec[1] += len(ivs)
            rec[2] += sum(d * (e - s + 1) for s, e, d in ivs)
            allint.extend((s, e) for s, e, _d in ivs)
        total_masked += merge(allint)
        buf.clear()

    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            try:
                div, query = float(f[1]), f[4]
                qbeg, qend = int(f[5]), int(f[6])
                entry = base_name(f[9])
            except ValueError:
                continue
            if query != cur_query:
                flush()
                cur_query = query
            buf.setdefault(entry, []).append((min(qbeg, qend), max(qbeg, qend), div))
    flush()
    return per_entry, total_masked


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("replicates", nargs="+", metavar="LABEL=FILE",
                    help="RepeatMasker .out file per read sample")
    ap.add_argument("--totals", required=True,
                    help="TSV with columns replicate, n_reads, total_bp")
    ap.add_argument("--groups", default=None,
                    help="groups.tsv from sat_group.py (adds family / superfamily)")
    ap.add_argument("--out", required=True, help="output per-entry abundance TSV")
    ap.add_argument("--summary", default=None, help="output per-replicate summary TSV")
    args = ap.parse_args()

    totals = {}
    with open(args.totals) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            totals[row["replicate"]] = (int(row["n_reads"]), int(row["total_bp"]))

    groups = {}
    if args.groups and os.path.isfile(args.groups):
        with open(args.groups) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                groups[row["seq_id"]] = row

    labels, per_rep, sat_frac = [], {}, {}
    for spec in args.replicates:
        if "=" not in spec:
            ap.error("replicate argument must look like LABEL=FILE, got %r" % spec)
        label, path = spec.split("=", 1)
        if label not in totals:
            ap.error("no total_bp for replicate %r in %s" % (label, args.totals))
        labels.append(label)
        entries, masked = scan(path)
        per_rep[label] = entries
        sat_frac[label] = 100.0 * masked / totals[label][1] if totals[label][1] else 0.0
        sys.stderr.write("[quantify] %s: %d entr(y/ies) hit, %.4f%% of %d bp masked\n"
                         % (label, len(entries), sat_frac[label], totals[label][1]))

    hit = {e for r in per_rep.values() for e in r}
    # Library entries are the variant representatives; the collapsed members are
    # only listed if something matched them (it should not, they are not masked
    # against), so the table stays one row per library entry.
    in_lib = {e for e, g in groups.items() if g.get("is_representative") == "1"}
    all_entries = sorted(hit | in_lib) if groups else sorted(hit)

    fields = (["seq_id", "monomer_len", "variant_group", "family_group",
               "superfamily_group"]
              + ["prop_pct_%s" % l for l in labels]
              + ["mean_prop_pct", "sd_prop_pct", "mean_divergence", "total_hits"])
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        for entry in all_entries:
            g = groups.get(entry, {})
            props, hits, div_num, div_den = [], 0, 0.0, 0
            for label in labels:
                rec = per_rep[label].get(entry)
                bp = rec[0] if rec else 0
                props.append(100.0 * bp / totals[label][1] if totals[label][1] else 0.0)
                if rec:
                    hits += rec[1]
                    div_num += rec[2]
                    div_den += rec[0]
            mean = sum(props) / len(props) if props else 0.0
            if len(props) > 1:
                var = sum((p - mean) ** 2 for p in props) / (len(props) - 1)
                sd = math.sqrt(var)
            else:
                sd = 0.0
            w.writerow([entry, g.get("monomer_len", "NA"), g.get("variant_group", "NA"),
                        g.get("family_group", "NA"), g.get("superfamily_group", "NA")]
                       + ["%.6f" % p for p in props]
                       + ["%.6f" % mean, "%.6f" % sd,
                          "%.2f" % (div_num / div_den) if div_den else "NA", hits])

    if args.summary:
        with open(args.summary, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["replicate", "n_reads", "total_bp", "satDNA_pct"])
            for label in labels:
                n, bp = totals[label]
                w.writerow([label, n, bp, "%.6f" % sat_frac[label]])

    sys.stderr.write("[quantify] %d library entr(y/ies) -> %s\n"
                     % (len(all_entries), args.out))


if __name__ == "__main__":
    main()
