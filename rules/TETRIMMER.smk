rule TETRIMMER:
    """Run TEtrimmer to clean and classify the RepeatModeler raw library against the genome."""
    input:
        rm_done = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}.repeatmodeler.done"),
        genome  = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna")
    output:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done")
    params:
        input_file         = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.fa"),
        output_dir         = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer"),
        pfam_dir           = config["PFAM_DIR"],
        max_cluster_number = 4,
        ext_step           = 2500,
        max_ext            = 10000,
        preset             = "conserved",
    log:
        os.path.join(LOG_DIR, "TEtrimmer", "tetrimmer_{species}.log")
    resources:
        cpus_per_task  = 20,
        mem_mb_per_cpu = 5000,
        runtime        = 720
    envmodules:
        "TEtrimmer/1.7.2-foss-2025a",
        "IQ-TREE/3.1.2-gompi-2025a"
    shell:
        """
        TEtrimmer \
            --input_file {params.input_file} \
            --genome_file {input.genome} \
            --output_dir {params.output_dir} \
            --pfam_dir {params.pfam_dir} \
            --num_threads {resources.cpus_per_task} \
            --preset {params.preset} \
            --max_cluster_number {params.max_cluster_number} \
            --ext_step {params.ext_step} \
            --max_ext {params.max_ext} \
            &> {log}
        touch {output.done}
        """
