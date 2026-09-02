import os
import glob
import gzip

GENOMES_DIR = config["GENOMES_DIR"]
GENOMES_DIR_DONE = config["GENOMES_DIR_DONE"]
GENOMES_SOURCE_DIR = config["GENOMES_SOURCE_DIR"]
LOG_DIR = config["LOG_DIR"]

# --- RepeatExplorer2 short-read branch ---
READS_DIR = config["READS_DIR"]
REPEATEXPLORER_SIF = config["REPEATEXPLORER_SIF"]
SEQTK = config["SEQTK"]
RE_TAXON = config["RE_TAXON"]
RE_COVERAGE = config["RE_COVERAGE"]
RE_COVERAGE_R2 = config["RE_COVERAGE_R2"]
RE_MINCL = config["RE_MINCL"]
RE_READLEN = config["RE_READLEN"]
RE_SEED = config["RE_SEED"]

SPECIES = [d for d in os.listdir(GENOMES_DIR_DONE)
            if os.path.isdir(os.path.join(GENOMES_DIR_DONE, d))]


def find_reads(species):
    """Find the paired-end read pair (R1, R2) for a species under READS_DIR.

    Tolerates both '*_1.fastq.gz'/'*_2.fastq.gz' and '*_1.fq.gz'/'*_2.fq.gz'
    naming. Returns (r1, r2) or None if no usable pair exists.
    """
    d = os.path.join(READS_DIR, species)
    for r1pat, r2pat in [("*_1.fastq.gz", "*_2.fastq.gz"), ("*_1.fq.gz", "*_2.fq.gz")]:
        r1 = sorted(glob.glob(os.path.join(d, r1pat)))
        r2 = sorted(glob.glob(os.path.join(d, r2pat)))
        if r1 and r2:
            return r1[0], r2[0]
    return None


def first_read_len(fq_gz):
    """Length of the first read's sequence in a gzipped FASTQ (fixed-length proxy)."""
    with gzip.open(fq_gz, "rt") as fh:
        fh.readline()                       # header
        return len(fh.readline().strip())   # sequence


# Species eligible for RepeatExplorer2: 
    # A genome dir
    # PE read pair
    # Reads with at least RE_READLEN long. 

SPECIES_WITH_READS = []
for _s in SPECIES:
    _rp = find_reads(_s)
    if _rp is None:
        continue                            # no PE pair -> silently skip
    if first_read_len(_rp[0]) < RE_READLEN:
        print(f"[RepeatExplorer] SKIP {_s}: reads shorter than {RE_READLEN} bp, "
              f"cannot enforce uniform length; not processed.")
        continue
    SPECIES_WITH_READS.append(_s)


# Species whose TEtrimmer consensus library already exists on disk.
# ONLY used by the `priority_all` convenience target, which is meant to pick up
# finished TEtrimmer runs without pulling TEtrimmer back into the DAG.
# Do NOT gate `rule all` on this: it is evaluated at parse time, so before the
# first TEtrimmer run it is empty and every downstream target silently
# disappears (that is what made `--until PRIORITY_TABLE` report 0 jobs).
# `rule all` uses SPECIES and lets the rule graph -- TETRIMMER -> REPEATMASKER
# -> RMSK_DIVERGENCE -> PRIORITY_TABLE, chained through the
# {species}.tetrimmer.done sentinel -- decide what still has to run.
SPECIES_WITH_TETRIMMER = [
    _s for _s in SPECIES
    if os.path.exists(os.path.join(GENOMES_DIR_DONE, _s, "TEtrimmer", "TEtrimmer_consensus_merged.fasta"))
]

include: "rules/1_RENAME_HEADERS.smk"
include: "rules/2_BUILD_RM_DATABASE.smk"
include: "rules/3_RUN_RM2.smk"
include: "rules/4_RECLASSIFY_RENAME.smk"
include: "rules/5_TETRIMMER.smk"
#include: "rules/6_REPEATEXPLORER.smk"
include: "rules/7_REPEATMASKER.smk"
include: "rules/8_PRIORITY_TABLE.smk"

rule all:
    input:
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),                            species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}.builddb.done"),                   species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}.repeatmodeler.done"),             species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}_rm1.0.fasta"),                    species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done"),            species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.divsum.html"),                    species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "Priority", "final_priority.table.tab"),             species=SPECIES)
#        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer.done"),  species=SPECIES_WITH_READS)
