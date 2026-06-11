import os

GENOMES_DIR = config["GENOMES_DIR"]
GENOMES_DIR_DONE = config["GENOMES_DIR_DONE"]
GENOMES_SOURCE_DIR = config["GENOMES_SOURCE_DIR"]
LOG_DIR = config["LOG_DIR"]
CURATEDLIB = config["CURATEDLIB"]

SPECIES = [d for d in os.listdir(GENOMES_DIR)
            if os.path.isdir(os.path.join(GENOMES_DIR, d))]

#include: "rules/rename_headers.smk"
#include: "rules/EDTA.smk"

include: "rules/RM_DATABASE.smk"
include: "rules/RM2.smk"
include: "rules/TETRIMMER.smk"
rule all:
    input:
        expand(os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}.builddb.done"), species=SPECIES),
        expand(os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}.repeatmodeler.done"), species=SPECIES),
        expand(os.path.join(GENOMES_DIR, "{species}", "TEtrimmer", "{species}.tetrimmer.done"), species=SPECIES)