#!/usr/bin/env python3
"""Group satDNA monomers into variants, families and superfamilies.

Applies the Ruiz-Ruano et al. (2016) classification criteria to the pairwise
table written by sat_homology.py, using single-linkage clustering at three
nested identity thresholds:

  >= --variant-id      (default 95%)  same variant  -> collapsed, one consensus kept
  >= --family-id       (default 80%)  different variants of the same family
  >= --superfamily-id  (default 50%)  different families of the same superfamily

The superfamily threshold is the soft one in the literature: the characid fish
satellitome uses >50%, the Chorthippus parallelus paper uses "any detectable
homology below 80%". Pass --superfamily-id 0 for the latter. Whichever you
use, state it in the methods.

Because the thresholds are nested and share one edge set, variants are always
contained in families and families in superfamilies.

Usage:
  sat_group.py --monomers sat.fasta --table sat.tsv --pairwise pairwise.tsv \
      --groups groups.tsv --nonredundant nonredundant.fasta
"""
import argparse
import csv
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


class UnionFind:
    """Minimal union-find over a fixed set of ids (single-linkage clustering)."""

    def __init__(self, items):
        self.parent = {i: i for i in items}

    def find(self, a):
        root = a
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[a] != root:          # path compression
            self.parent[a], a = root, self.parent[a]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

    def groups(self):
        out = {}
        for i in self.parent:
            out.setdefault(self.find(i), []).append(i)
        return out


def cluster(ids, edges, threshold):
    """{seq_id: group_label} for single-linkage clustering at `threshold`.

    Group labels are numbered by the alphabetically first member so that the
    same input always yields the same labels.
    """
    uf = UnionFind(ids)
    for (a, b), ident in edges.items():
        if ident >= threshold:
            uf.union(a, b)
    members = sorted(uf.groups().values(), key=lambda m: sorted(m)[0])
    return {sid: n for n, m in enumerate(members, 1) for sid in m}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--monomers", required=True)
    ap.add_argument("--table", required=True, help="metadata TSV from collect_tarean_satellites.py")
    ap.add_argument("--pairwise", required=True, help="TSV from sat_homology.py")
    ap.add_argument("--groups", required=True, help="output group assignment TSV")
    ap.add_argument("--nonredundant", required=True,
                    help="output FASTA with one consensus per variant")
    ap.add_argument("--variant-id", type=float, default=95.0)
    ap.add_argument("--family-id", type=float, default=80.0)
    ap.add_argument("--superfamily-id", type=float, default=50.0,
                    help="0 = any detectable RepeatMasker homology")
    ap.add_argument("--min-cov", type=float, default=0.5,
                    help="minimum fraction of the shorter monomer that must align "
                         "for a hit to count as homology (default 0.5)")
    args = ap.parse_args()

    seqs = {h.split()[0]: s for h, s in read_fasta(args.monomers)}
    meta = {}
    with open(args.table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["seq_id"]] = row

    # Symmetrise the directional pairwise table: keep the better of the two runs.
    edges = {}
    with open(args.pairwise) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                ident, cov = float(row["identity"]), float(row["coverage"])
            except (KeyError, ValueError):
                continue
            if cov < args.min_cov:
                continue
            key = tuple(sorted((row["query"], row["subject"])))
            if key[0] == key[1] or key[0] not in seqs or key[1] not in seqs:
                continue
            edges[key] = max(edges.get(key, 0.0), ident)

    ids = sorted(seqs)
    variant = cluster(ids, edges, args.variant_id)
    family = cluster(ids, edges, args.family_id)
    superfam = cluster(ids, edges, args.superfamily_id if args.superfamily_id > 0
                       else 0.0)

    def rank_key(sid):
        """Representative preference: TAREAN rank, then cluster size, then length."""
        m = meta.get(sid, {})
        try:
            tarean = int(m.get("tarean_rank", "9"))
        except ValueError:
            tarean = 9
        try:
            size = float(m.get("re_cluster_size", "0") or 0)
        except ValueError:
            size = 0.0
        return (tarean, -size, -len(seqs[sid]), sid)

    reps = {}
    for sid, v in variant.items():
        if v not in reps or rank_key(sid) < rank_key(reps[v]):
            reps[v] = sid

    fields = ["seq_id", "round", "cluster", "tarean_rank", "confidence", "monomer_len",
              "re_cluster_size", "re_cluster_prop_pct", "variant_group", "family_group",
              "superfamily_group", "is_representative", "representative"]
    with open(args.groups, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for sid in ids:
            m = meta.get(sid, {})
            w.writerow({
                "seq_id": sid,
                "round": m.get("round", "NA"),
                "cluster": m.get("cluster", "NA"),
                "tarean_rank": m.get("tarean_rank", "NA"),
                "confidence": m.get("confidence", "NA"),
                "monomer_len": len(seqs[sid]),
                "re_cluster_size": m.get("re_cluster_size", "NA"),
                "re_cluster_prop_pct": m.get("re_cluster_prop_pct", "NA"),
                "variant_group": "V%03d" % variant[sid],
                "family_group": "F%03d" % family[sid],
                "superfamily_group": "SF%03d" % superfam[sid],
                "is_representative": int(reps[variant[sid]] == sid),
                "representative": reps[variant[sid]],
            })

    with open(args.nonredundant, "w") as fh:
        for v in sorted(reps):
            sid = reps[v]
            fh.write(">%s variant=V%03d family=F%03d superfamily=SF%03d len=%d\n%s\n"
                     % (sid, v, family[sid], superfam[sid], len(seqs[sid]), seqs[sid]))

    sys.stderr.write(
        "[sat_group] %d monomer(s) -> %d variant(s) (>=%.0f%%), %d famil(y/ies) "
        "(>=%.0f%%), %d superfamil(y/ies) (%s)\n"
        % (len(ids), len(set(variant.values())), args.variant_id,
           len(set(family.values())), args.family_id, len(set(superfam.values())),
           "any homology" if args.superfamily_id <= 0
           else ">=%.0f%%" % args.superfamily_id))


if __name__ == "__main__":
    main()
