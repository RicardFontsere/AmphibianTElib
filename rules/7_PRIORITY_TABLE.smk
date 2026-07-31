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
        pfamdb   = config["PFAM_DIR"],
        repbase  = config["REPBASE_DB"],
        pfamscan = config["PFAMSCAN_DIR"],
        script   = os.path.join(workflow.basedir, "scripts", "priority_table.sh"),
    log:
        os.path.join(LOG_DIR, "PriorityTable", "priority_{species}.log")
    resources:
        cpus_per_task  = 2,
        mem_mb_per_cpu = 1000,
        runtime        = 20
    envmodules:
        "CD-HIT/4.8.1-GCC-13.3.0",
        "BLAST+/2.16.0-gompi-2024a",
        "EMBOSS/6.6.0-foss-2024a",
        "BioPerl/1.7.8-GCCcore-13.3.0",
        "HMMER/3.4-gompi-2024a"
    shell:
        """
        mkdir -p {params.workdir} $(dirname {log})
        cd {params.workdir}
        bash {params.script} \
            {input.library} \
            {input.genome} \
            {params.pfamdb} \
            {params.repbase} \
            {resources.cpus_per_task} \
            {params.pfamscan} \
            &> {log}
        """
