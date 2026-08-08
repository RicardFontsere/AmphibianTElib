rule rename_fasta_headers:
    """Truncate FASTA headers at the first dot for EDTA compatibility.
       >NC_090921.1  →  >NC_090921
    """
    output:
        fasta = os.path.join(GENOMES_DIR, "{species}", "{species}_headers.fna")
    params:
        species_dir = os.path.join(GENOMES_DIR, "{species}")
    log:
        os.path.join(LOG_DIR, "rename_fasta_headers", "rename_fasta_headers_{species}.log")
    shell:
        """
        GENOME=$(find {params.species_dir} -name '*.fna' ! -name '*cds*' | head -1) #Earlier runs contained cds files, not anymore could be cleaned up.
        awk '/^>/ {{ sub(/\\..*/, ""); }} {{ print }}' "$GENOME" > {output.fasta}
        echo "rename_fasta_headers completed for {wildcards.species}: $GENOME -> {output.fasta}" > {log}
        """
 