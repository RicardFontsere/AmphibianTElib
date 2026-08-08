rule priority_all:
    """Convenience target: build ONLY the priority tables for species that already
    have a TEtrimmer library. Run `snakemake priority_all ...` to avoid pulling the
    rest of `rule all` (which would try to build TEtrimmer etc. for every species)."""
    input:
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority", "final_priority.table.tab"),
               species=SPECIES_WITH_TETRIMMER)


rule PRIORITY_TABLE:
    """Post-process the TEtrimmer library into a per-family priority table:
    cd-hit-est dereplication, consensus length, genome copy-number, RepBase
    novelty screen (against a pre-built RepBase blast db) and Pfam domains,
    merged into final_priority.table.tab. Runs in place in TEtrimmer/priority/."""
    input:
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "TEtrimmer_consensus_merged.fasta"),
        genome  = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
    output:
        table = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority", "final_priority.table.tab"),
    params:
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority"),
        cdddir  = config["CDD_DIR"],
        repbase = config["REPBASE_DB"],
        script  = os.path.join(workflow.basedir, "scripts", "priority_table.sh"),
    log:
        os.path.join(LOG_DIR, "PriorityTable", "priority_{species}.log")
    resources:
        cpus_per_task  = 4,
        mem_mb_per_cpu = 10000,
        runtime        = 200
    envmodules:
        "CD-HIT/4.8.1-GCC-13.3.0",
        "BLAST+/2.16.0-gompi-2024a",
        "GCC/13.3.0"
    shell:
        """
        mkdir -p {params.workdir} $(dirname {log})
        cd {params.workdir}
        bash {params.script} \
            {input.library} \
            {input.genome} \
            {params.cdddir} \
            {params.repbase} \
            {resources.cpus_per_task} \
            &> {log}
        """
