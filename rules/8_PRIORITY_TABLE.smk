rule priority_all:
    """Convenience target: build ONLY the priority tables for species that already
    have a TEtrimmer library. Run `snakemake priority_all ...` to avoid pulling the
    rest of `rule all` (which would try to build TEtrimmer etc. for every species)."""
    input:
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "Priority", "final_priority.table.tab"),
               species=SPECIES_WITH_TETRIMMER)



rule PRIORITY_TABLE:
    """Post-process the TEtrimmer library into a per-family priority table:
    cd-hit-est dereplication, consensus length, genome occupancy and Kimura age
    from the RepeatMasker .divsum, defragmented/full-length copy counts from the
    RepeatMasker .out, RepBase novelty screen, CDD domains and a TRF satellite
    screen, merged into final_priority.table.tab. Runs in place in Priority/,
    at the same depth as RMDB/ and RMSK/."""
    input:
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "TEtrimmer_consensus_merged.fasta"),
        divsum  = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.align.divsum"),
        rmout   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.out"),
    output:
        table = os.path.join(GENOMES_DIR_DONE, "{species}", "Priority", "final_priority.table.tab"),
    params:
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "Priority"),
        cdddir  = config["CDD_DIR"],
        repbase = config["REPBASE_DB"],
        script  = os.path.join(workflow.basedir, "scripts", "priority_table_withTRF.sh"),
    log:
        os.path.join(LOG_DIR, "PriorityTable", "priority_{species}.log")
    resources:
        cpus_per_task  = 10,
        mem_mb_per_cpu = 8000,
        runtime        = 600
    envmodules:
        "CD-HIT/4.8.1-GCC-13.3.0",
        "BLAST+/2.16.0-gompi-2024a",
        "GCC/13.3.0",
        "TRF/4.09.1-GCCcore-13.3.0"
    shell:
        """
        mkdir -p {params.workdir} $(dirname {log})
        cd {params.workdir}
        bash {params.script} \
            {input.library} \
            {input.divsum} \
            {input.rmout} \
            {params.cdddir} \
            {params.repbase} \
            {resources.cpus_per_task} \
            &> {log}
        """
