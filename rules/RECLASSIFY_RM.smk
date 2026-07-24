rule RECLASSIFY_RM:
    """Re-run RepeatClassifier on the RepeatModeler library using a complete
    RepeatMasker/Dfam library. The original de-novo run classified against an
    incomplete RepeatMasker.lib (missing Dfam partitions), leaving many families
    Unknown/misclassified. This linearises the raw families FASTA and reclassifies
    it in place -- it does NOT re-run de-novo discovery."""
    input:
        rm_done  = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}.repeatmodeler.done"),
        families = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}-families.fa"),
        stk      = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}-families.stk")
    output:
        reclassified = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}-families.reclassified.fa")
    params:
        workdir = os.path.join(GENOMES_DIR, "{species}", "RMDB", "reclassify"),
        libdir  = config.get("RM_LIBDIR", "")   # optional: complete RepeatMasker Libraries dir
    log:
        os.path.join(LOG_DIR, "Reclassify", "reclassify_{species}.log")
    resources:
        cpus_per_task  = 16,
        mem_mb_per_cpu = 4000,
        runtime        = 2880
    envmodules:
        "RepeatModeler/2.0.8-foss-2025a",
        "Perl-bundle-CPAN/5.40.0-GCCcore-14.2.0",
        "BioPerl/1.7.8-GCCcore-14.2.0"
    shell:
        """
        mkdir -p {params.workdir} $(dirname {log})

        # Optional: point RepeatClassifier at a complete RepeatMasker/Dfam library.
        # If RM_LIBDIR is empty, the loaded module's default library is used.
        if [ -n "{params.libdir}" ]; then export LIBDIR="{params.libdir}"; fi

        # 1. Linearise the RepeatModeler library to single-line FASTA
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
        """
