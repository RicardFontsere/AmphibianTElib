rule TETRIMMER:
    """Run TEtrimmer on raw library against the genome."""
    input:
        genome  = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}_rm1.0.fasta"),
    output:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done"),
    params:
        input_file         = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.fa"),
        output_dir         = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer"),
        pfam_dir           = config["PFAM_DIR"],
        preset             = "conserved",
    log:
        os.path.join(LOG_DIR, "TEtrimmer", "tetrimmer_{species}.log")
    resources:
        cpus_per_task  = 16,
        mem_mb_per_cpu = 28000,
        runtime        = 2200
    envmodules:
        "TEtrimmer/1.7.2-foss-2025a",
    shell:
        """
        TEtrimmer \
            --input_file {input.library} \
            --genome_file {input.genome} \
            --output_dir {params.output_dir} \
            --pfam_dir {params.pfam_dir} \
            --num_threads {resources.cpus_per_task} \
            --preset {params.preset} \
            &> {log}
        touch {output.done}
        """
