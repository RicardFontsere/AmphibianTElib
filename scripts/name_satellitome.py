#!/usr/bin/env python3
"""Name the satDNA families and write the final satellitome library.

Naming follows Ruiz-Ruano et al. (2016): species abbreviation + "Sat" + a
catalogue number assigned in order of decreasing abundance + the consensus
monomer length, e.g. LmiSat06-185. Families made of several variants (80-95%
identity to each other) get a letter suffix per variant, e.g. RteSat03-176A
and RteSat03-180B; families with a single variant get no letter.

Families are ranked by the summed abundance of their variants, as re-estimated
by quantify_satellitome.py -- not by RepeatExplorer cluster size.

Usage:
  name_satellitome.py --species Rana_temporaria --nonredundant nr.fasta \
      --groups groups.tsv --abundance abundance.tsv \
      --out-fasta lib.fasta --catalogue catalogue.tsv
"""
import argparse
import csv
import string
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


def abbreviate(species):
    """'Rana_temporaria' -> 'Rte' (genus initial + first two letters of epithet)."""
    parts = species.replace(" ", "_").split("_")
    genus = parts[0] if parts else species
    epithet = parts[1] if len(parts) > 1 else ""
    return (genus[:1].upper() + epithet[:2].lower()) or species[:3]


def letters(n):
    """0 -> A, 25 -> Z, 26 -> AA, ..."""
    out = ""
    n += 1
    while n:
        n, r = divmod(n - 1, 26)
        out = string.ascii_uppercase[r] + out
    return out


def at_pct(seq):
    acgt = sum(seq.count(b) for b in "ACGT")
    return 100.0 * (seq.count("A") + seq.count("T")) / acgt if acgt else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True, help="Genus_epithet")
    ap.add_argument("--abbrev", default="auto", help="library prefix (default: derived)")
    ap.add_argument("--nonredundant", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--abundance", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--catalogue", required=True)
    ap.add_argument("--version", default="1.0")
    ap.add_argument("--min-abundance", type=float, default=0.0,
                    help="drop variants whose mean abundance is below this %% "
                         "(default 0.0 = keep everything TAREAN called)")
    args = ap.parse_args()

    prefix = abbreviate(args.species) if args.abbrev == "auto" else args.abbrev

    seqs = {h.split()[0]: s.upper() for h, s in read_fasta(args.nonredundant)}

    rows, members = {}, {}
    with open(args.groups) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            members.setdefault(row["representative"], []).append(row["seq_id"])
            if row.get("is_representative") == "1":
                rows[row["seq_id"]] = row

    abundance = {}
    with open(args.abundance) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            abundance[row["seq_id"]] = row

    def prop(sid):
        try:
            return float(abundance.get(sid, {}).get("mean_prop_pct", 0.0))
        except ValueError:
            return 0.0

    entries = [sid for sid in seqs if sid in rows and prop(sid) >= args.min_abundance]
    dropped = len(seqs) - len(entries)

    # Rank families by total abundance, then variants within family by abundance.
    by_family = {}
    for sid in entries:
        by_family.setdefault(rows[sid]["family_group"], []).append(sid)
    order = sorted(by_family, key=lambda f: (-sum(prop(s) for s in by_family[f]), f))

    catalogue = []
    for n, fam in enumerate(order, 1):
        variants = sorted(by_family[fam], key=lambda s: (-prop(s), s))
        multi = len(variants) > 1
        for i, sid in enumerate(variants):
            seq = seqs[sid]
            name = "%sSat%02d-%d%s" % (prefix, n, len(seq),
                                       letters(i) if multi else "")
            ab = abundance.get(sid, {})
            catalogue.append({
                "name": name,
                "family_no": n,
                "variant": letters(i) if multi else "",
                "superfamily": rows[sid]["superfamily_group"],
                "monomer_len": len(seq),
                "at_pct": "%.1f" % at_pct(seq),
                "abundance_pct": "%.6f" % prop(sid),
                "sd_pct": ab.get("sd_prop_pct", "NA"),
                "divergence": ab.get("mean_divergence", "NA"),
                "tarean_confidence": rows[sid].get("confidence", "NA"),
                "re_round": rows[sid].get("round", "NA"),
                "re_cluster": rows[sid].get("cluster", "NA"),
                "source_id": sid,
                "collapsed_ids": ",".join(sorted(members.get(sid, [sid]))),
                "seq": seq,
            })

    with open(args.out_fasta, "w") as fh:
        for rec in catalogue:
            fh.write(">%s#Satellite superfamily=%s abundance=%s%% monomer=%d source=%s\n%s\n"
                     % (rec["name"], rec["superfamily"], rec["abundance_pct"],
                        rec["monomer_len"], rec["source_id"], rec["seq"]))

    fields = ["name", "family_no", "variant", "superfamily", "monomer_len", "at_pct",
              "abundance_pct", "sd_pct", "divergence", "tarean_confidence", "re_round",
              "re_cluster", "source_id", "collapsed_ids"]
    with open(args.catalogue, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(catalogue)

    total = sum(float(r["abundance_pct"]) for r in catalogue)
    sys.stderr.write("[name_satellitome] %s v%s: %d famil(y/ies), %d librar(y/ies) "
                     "entr(y/ies), %.4f%% of the genome%s\n"
                     % (prefix, args.version, len(order), len(catalogue), total,
                        ", %d below --min-abundance" % dropped if dropped else ""))


if __name__ == "__main__":
    main()
