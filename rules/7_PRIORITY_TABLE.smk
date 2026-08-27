rule priority_all:
    """Convenience target: build ONLY the priority tables (and the linearised
    lib2.0 FASTA) for species that already have a TEtrimmer library. Run
    `snakemake priority_all ...` to avoid pulling the rest of `rule all`
    (which would try to build TEtrimmer etc. for every species)."""
    input:
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority", "final_priority.table.tab"),
               species=SPECIES_WITH_TETRIMMER),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}_lib2.0.fa"),
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
        cpus_per_task  = 10,
        mem_mb_per_cpu = 8000,
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


rule TETRIMMER_LIB2:
    """Linearise the TEtrimmer consensus library into the final `lib2.0` FASTA:
    every record is folded onto a single sequence line and each header is
    truncated at the first whitespace, so downstream tools see clean IDs.
    Output is named after the species (genus_epithet), e.g.
    Bufotes_viridis -> Bufotes_viridis_lib2.0.fa."""
    input:
        done    = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done"),
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "TEtrimmer_consensus_merged.fasta"),
    output:
        lib = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}_lib2.0.fa"),
    log:
        os.path.join(LOG_DIR, "PriorityTable", "lib2.0_{species}.log")
    resources:
        cpus_per_task  = 1,
        mem_mb_per_cpu = 4000,
        runtime        = 30
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.lib})
        awk '/^>/ {{sub(/[[:space:]].*/, ""); printf("\\n%s\\n", $0); next}} {{printf("%s", $0)}} END {{printf("\\n")}}' \
            < {input.library} \
          | awk 'NR > 1' > {output.lib} 2> {log}
        """
