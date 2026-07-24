rule edta:
    """Run EDTA on the renamed-headers genome with CDS and curated library."""
    input:
        genome = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
        rmlib  = os.path.join(GENOMES_DIR_DONE, "{species}", "RMDB", "{species}-families.fa")
    output:
        helitron = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_helitron.denovo")
    params:
        species_dir = os.path.join(GENOMES_DIR_DONE, "{species}"),
    log:
        os.path.join(LOG_DIR, "EDTA_Helitron", "EDTA_Helitron_{species}.log")
    resources:
        cpus_per_task=16,
        mem_mb_per_cpu=3000,
        runtime=2888
    shell:
        """
        cd {params.species_dir}
        PYTHONNOUSERSITE=1 apptainer exec /apps/brussel/containers/edta/EDTA-2.3.0--hdfd78af_0.sif bash -c \
            'export LC_ALL=C LANG=C LANGUAGE=C && unset -f which && EDTA_raw.pl \
            --genome {input.genome} \
            --species others \
            --type helitron \
            --rmlib {input.rmlib}
            --threads {resources.cpus_per_task}' 2> {log}
        touch {output.helitron}
        """
