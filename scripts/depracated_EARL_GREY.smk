rule EARL_GREY:
    """Run EarlGrey TE annotation on the renamed-headers genome with a curated library."""
    input:
        genome = os.path.join(GENOMES_DIR, "{species}", "{species}_headers.fna")
    output:
        log_file = os.path.join(GENOMES_DIR, "{species}", "EarlGrey", "{species}EarlGrey.log")
    params:
        species_dir    = os.path.join(GENOMES_DIR, "{species}"),
        genomes_source = os.path.join(GENOMES_SOURCE_DIR, "{species}"),
        te_dir         = os.path.dirname(os.path.dirname(CURATEDLIB)),
        curatedlib     = CURATEDLIB,
        flank          = 3000,
        iters          = 3
    log:
        os.path.join(LOG_DIR, "EarlGrey", "earlgrey_{species}.log")
    resources:
        cpus_per_task = 32,
        mem_mb_per_cpu = 3000,
        runtime = 3000
    shell:
        """
        apptainer exec \
            --bind {params.species_dir}:/data/ \
            --bind {params.genomes_source}:/genomes \
            --bind {params.te_dir}:{params.te_dir} \
            /apps/brussel/containers/earlgrey/EarlGrey-7.2.1.sif sh -c \
            "unset -f which && earlGrey \
                -g /genomes/{wildcards.species}_headers.fna \
                -s {wildcards.species} \
                -o /data/ \
                -t {resources.cpus_per_task} \
                -i {params.iters} \
                -l {params.curatedlib} \
                -f {params.flank} " &> {log}
        """

