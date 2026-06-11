rule TETRIMMER:
    """Run TEtrimmer to clean and classify the RepeatModeler raw library against the genome."""
    input:
        rm_done = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}.repeatmodeler.done"),
        genome  = os.path.join(GENOMES_DIR, "{species}", "{species}_headers.fna")
    output:
        done = os.path.join(GENOMES_DIR, "{species}", "TEtrimmer", "{species}.tetrimmer.done")
    params:
        input_file        = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}-families.fa"),
        output_dir        = os.path.join(GENOMES_DIR, "{species}", "TEtrimmer"),
        pfam_dir          = config["PFAM_DIR"],
        max_msa_lines     = 100,
        top_msa_lines     = 90,
        max_cluster_number = 3,
        ext_thr           = 10000,
        ext_check_win     = 2500
    log:
        os.path.join(LOG_DIR, "TEtrimmer", "tetrimmer_{species}.log")
    resources:
        cpus_per_task  = 20,
        mem_mb_per_cpu = 5000,
        runtime        = 720
    envmodules:
        "TEtrimmer/1.5.4-foss-2024a"
    shell:
        """
        TEtrimmer \
            --input_file {params.input_file} \
            --genome_file {input.genome} \
            --output_dir {params.output_dir} \
            --pfam_dir {params.pfam_dir} \
            --max_msa_lines {params.max_msa_lines} \
            --top_msa_lines {params.top_msa_lines} \
            --num_threads {resources.cpus_per_task} \
            --max_cluster_number {params.max_cluster_number} \
            --ext_thr {params.ext_thr} \
            --ext_check_win {params.ext_check_win} \
            --classify_all \
            &> {log}
        touch {output.done}
        """
