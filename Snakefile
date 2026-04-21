import os

configfile: "config.yaml"

GENOMES_DIR = config["genomes_dir"]
LOG_DIR = config["log_dir"]
CURATEDLIB = config["curatedlib"]
SPECIES = [d for d in os.listdir(GENOMES_DIR)
            if os.path.isdir(os.path.join(GENOMES_DIR, d))]

include: "rules/rename_headers.smk"
include: "rules/EDTA.smk"
#include: "rules/EDTA_TIR.smk"

rule all:
    input:
        expand(os.path.join(GENOMES_DIR, "{species}", "{species}.denovo"), species=SPECIES)

