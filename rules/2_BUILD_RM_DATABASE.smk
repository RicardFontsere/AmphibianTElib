rule REPEATMODELER_DB:
    """Build a RepeatModeler BLAST database from the renamed-headers genome."""
    input:
        genome = os.path.join(GENOMES_DIR, "{species}", "{species}_headers.fna")
    output:
        marker = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}.builddb.done")
    params:
        rmdb_dir  = os.path.join(GENOMES_DIR, "{species}", "RMDB"),
        db_prefix = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}")
    log:
        os.path.join(LOG_DIR, "RepeatModeler", "repeatmodeler_db_{species}.log")
    resources:
        cpus_per_task  = 4,
        mem_mb_per_cpu = 2000,
        runtime        = 120
    envmodules:
        "RepeatModeler/2.0.8-foss-2025a",
        "Perl-bundle-CPAN/5.40.0-GCCcore-14.2.0",
        "BioPerl/1.7.8-GCCcore-14.2.0"
    shell:
        """
        mkdir -p {params.rmdb_dir}
        BuildDatabase -name {params.db_prefix} {input.genome} &> {log}
        touch {output.marker}
        """