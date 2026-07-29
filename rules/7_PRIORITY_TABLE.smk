rule PRIORITY_TABLE:
    """Post-process the TEtrimmer library into a per-family priority table:
    cd-hit-est dereplication, consensus length, genome copy-number, RepBase
    novelty screen (against a pre-built RepBase blast db) and Pfam domains,
    merged into final_priority.table.tab. Runs in place in TEtrimmer/priority/."""
    input:
        tetrimmer_done = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done"),
        genome         = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
    output:
        table = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority", "final_priority.table.tab"),
    params:
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "TEtrimmer_consensus_merged.fasta"),
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority"),
        pfamdb  = config["PFAM_DIR"],
        repbase = config["REPBASE_DB"],
        script  = os.path.join(workflow.basedir, "scripts", "priority_table.sh"),
    log:
        os.path.join(LOG_DIR, "PriorityTable", "priority_{species}.log")
    resources:
        cpus_per_task  = 8,
        mem_mb_per_cpu = 4000,
        runtime        = 720
    envmodules:
        # TODO: confirm/replace these module strings with the ones on your cluster
        "CD-HIT/4.8.1-GCC-13.3.0",
        "BLAST+/2.16.0-gompi-2024a",
        "EMBOSS/6.6.0-foss-2024a",
        "PfamScan/1.6-foss-2024a",
    shell:
        """
        mkdir -p {params.workdir} $(dirname {log})
        cd {params.workdir}
        bash {params.script} \
            {params.library} \
            {input.genome} \
            {params.pfamdb} \
            {params.repbase} \
            {resources.cpus_per_task} \
            &> {log}
        """
