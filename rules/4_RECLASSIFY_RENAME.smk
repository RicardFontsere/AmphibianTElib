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
    (folds in the former RENAME_RM_LIB rule)."""
    input:
        families = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.fa"),
        stk      = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.stk")
    output:
        reclassified = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.reclassified.fa"),
        renamed      = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}_rm1.0.fasta")
    params:
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "reclassify"),
        prefix  = species_prefix,
        script  = os.path.join(workflow.basedir, "scripts", "renameRMDLconsensi.pl"),
        tmp     = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}_rm1.0_temp.fasta")
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
        mkdir -p {params.workdir} $(dirname {log})
        #RepeatModeler library to single-line FASTA
        awk '/^>/ {{printf("\\n%s\\n",$0);next;}} {{printf("%s",$0);}} END {{printf("\\n");}}' \
            < {input.families} \
            | awk 'NR > 1' > {params.workdir}/consensi.fa 2> {log}

        # 2. Bring the Stockholm seed alignments alongside (IDs unchanged by linearisation)
        cp {input.stk} {params.workdir}/families.stk

        # 3. Reclassify -- RepeatClassifier writes consensi.fa.classified next to the input
        cd {params.workdir} && \
        RepeatClassifier \
            -consensi consensi.fa \
            -stockholm families.stk \
            -threads {resources.cpus_per_task} \
            >> {log} 2>&1

        # 4. Publish the reclassified single-line library
        cp {params.workdir}/consensi.fa.classified {output.reclassified}

        # 5. Rename headers with the short species prefix 
        awk '/^>/ {{printf("\\n%s\\n",$0);next;}} {{printf("%s",$0);}} END {{printf("\\n");}}' \
            < {output.reclassified} \
            | awk 'NR > 1' > {params.tmp} 2>> {log}
        perl {params.script} {params.tmp} {params.prefix} {output.renamed} 2>> {log}
        rm -f {params.tmp}
        echo "Finished renaming RepeatModeler library for {wildcards.species} " >> {log}
        """
