rule edta_tir:
    """Run EDTA on the renamed-headers genome with CDS and curated library."""
    input:
        genome = os.path.join(GENOMES_DIR, "{species}", "{species}_headers.fna")
    output:
        denovo = os.path.join(GENOMES_DIR, "{species}", "{species}.denovo")
    params:
        species_dir = os.path.join(GENOMES_DIR, "{species}"),
        curatedlib = CURATEDLIB
    log:
        os.path.join(LOG_DIR, "edta", "edta_{species}.log")
    resources:
        cpus_per_task=16,
        mem_mb_per_cpu=25000,
        runtime=7200
    shell:
        """
        CDS=$(find {params.species_dir} -name '*cds*.fna' | head -1)
        cd {params.species_dir}
        PYTHONNOUSERSITE=1 apptainer exec /apps/brussel/containers/edta/EDTA-2.3.0--hdfd78af_0.sif bash -c \
            'export LC_ALL=C LANG=C LANGUAGE=C && unset -f which && EDTA.pl \
            --genome {input.genome} \
            --cds '"$CDS"' \
            --curatedlib {params.curatedlib} \
            --overwrite 0 \
            --sensitive 0 \
            --anno 0 \
            --evaluate 0 \
            --threads {resources.cpus_per_task} > {output.denovo}' 2> {log}
        """
