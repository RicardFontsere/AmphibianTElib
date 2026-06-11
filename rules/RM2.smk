rule REPEATMODELER:
    """Run RepeatModeler2 with LTR structural discovery on the prebuilt database."""
    input:
        db_done = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}.builddb.done")
    output:
        done = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}.repeatmodeler.done")
    params:
        rmdb_dir  = os.path.join(GENOMES_DIR, "{species}", "RMDB"),
        db_prefix = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}")
    log:
        os.path.join(LOG_DIR, "RepeatModeler", "repeatmodeler_{species}.log")
    resources:
        cpus_per_task  = 64,
        mem_mb_per_cpu = 6000,
        runtime        = 7200
    envmodules:
        "RepeatModeler/2.0.8-foss-2025a",
        "Perl-bundle-CPAN/5.40.0-GCCcore-14.2.0",
        "BioPerl/1.7.8-GCCcore-14.2.0"
    shell:
        """
        cd {params.rmdb_dir} && \
        RepeatModeler \
            -database {wildcards.species} \
            -threads {resources.cpus_per_task} \
            -LTRStruct \
            -srand 1337 \
            -genomeSampleSizeMax 810000000 \
            &> {log}
        touch {output.done}
        """