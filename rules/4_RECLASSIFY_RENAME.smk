def species_prefix(wildcards):
    """Build a short library prefix from a 'genus_epithet' species name.

    First 3 letters of the genus (lower-case) + first 3 letters of the
    epithet with its initial capitalised, e.g. 'Ranitomeya_variabilis' -> 'ranVar'.
    """
    genus, epithet = wildcards.species.split("_", 1)
    return genus[:3].lower() + epithet[:3].capitalize()


rule RECLASSIFY_RM:
    """Re-run RepeatClassifier (temporary fix as previour RepeatMasker lacked libraries),
    then rename the reclassified library headers with a short species prefix
    (folds in the former RENAME_RM_LIB rule). Everything runs in place inside the
    species RMDB/ directory -- no copy-in / copy-out."""
    input:
        families = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.fa"),
        stk      = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.stk")
    output:
        reclassified = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.reclassified.fa"),
        renamed      = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}_rm1.0.fasta")
    params:
        rmdb_dir = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB"),
        prefix   = species_prefix,
        script   = os.path.join(workflow.basedir, "scripts", "renameRMDLconsensi.pl")
    log:
        os.path.join(LOG_DIR, "Reclassify", "reclassify_{species}.log")
    resources:
        cpus_per_task  = 8,
        mem_mb_per_cpu = 1000,
        runtime        = 1000
    envmodules:
        "RepeatModeler/2.0.8-foss-2025a",
        "Perl-bundle-CPAN/5.40.0-GCCcore-14.2.0",
        "BioPerl/1.7.8-GCCcore-14.2.0"
    shell:
        """
        mkdir -p $(dirname {log})
        cd {params.rmdb_dir}

        # 1. Linearise the RepeatModeler library to single-line FASTA (written in place)
        awk '/^>/ {{printf("\\n%s\\n",$0);next;}} {{printf("%s",$0);}} END {{printf("\\n");}}' \
            < {input.families} \
            | awk 'NR > 1' > {output.reclassified} 2> {log}

        # 2. Reclassify in place -- RepeatClassifier writes <consensi>.classified next to it.
        #    Stockholm alignments are read directly from RMDB/ (IDs unchanged by linearisation).
        RepeatClassifier \
            -consensi {output.reclassified} \
            -stockholm {input.stk} \
            -threads {resources.cpus_per_task} \
            >> {log} 2>&1

        # 3. RepeatClassifier appended .classified -> keep it as the reclassified library
        mv {output.reclassified}.classified {output.reclassified}

        # 4. Rename headers with the short species prefix (formerly RENAME_RM_LIB).
        #    The reclassified library is already single-line, so no re-linearisation needed.
        perl {params.script} {output.reclassified} {params.prefix} {output.renamed} 2>> {log}
        echo "Finished reclassifying + renaming RepeatModeler library for {wildcards.species} " >> {log}
        """
