import os

GENOMES_DIR = config["GENOMES_DIR"]
GENOMES_DIR_DONE = config["GENOMES_DIR_DONE"]
GENOMES_SOURCE_DIR = config["GENOMES_SOURCE_DIR"]
LOG_DIR = config["LOG_DIR"]
CURATEDLIB = config["CURATEDLIB"]

SPECIES = [d for d in os.listdir(GENOMES_DIR_DONE)
            if os.path.isdir(os.path.join(GENOMES_DIR_DONE, d))]

#include: "rules/rename_headers.smk"
#include: "rules/EDTA_HELITRON.smk"
include: "rules/RM_DATABASE.smk"
include: "rules/RM2.smk"
include: "rules/RENAME_RM.smk"
include: "rules/TETRIMMER.smk"
rule all:
    input:
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}.builddb.done"), species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}.repeatmodeler.done"), species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}_rm1.0.fasta"), species=SPECIES),
        #expand(os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_helitron.denovo"), species=SPECIES),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done"), species=SPECIES)
