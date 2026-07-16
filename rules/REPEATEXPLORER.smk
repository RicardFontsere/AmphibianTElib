def repeatexplorer_reads(wildcards):
    """Return the paired-end read pair (R1, R2) for a species.

    Species without a usable pair are filtered out of SPECIES_WITH_READS before
    the DAG is built, so this always finds a pair for the species it is called on.
    """
    r1, r2 = find_reads(wildcards.species)
    return {"r1": r1, "r2": r2}


rule SHORT_READ_PREP:
    """Prepare short reads for RepeatExplorer2: fastp QC, coverage-based subsample,
    uniform hard-trim, coded FASTA rename and interleave of the mates.

    Genome size (for the coverage target) is summed from the genome .fai index.
    Everything runs inside the species RepeatExplorer/ directory, and every
    intermediate file is prefixed with the species name to avoid confusion.
    """
    input:
        unpack(repeatexplorer_reads),
        genome = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
    output:
        interleaved = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_reads_interleaved.fasta"),
    params:
        seqtk    = SEQTK,
        coverage = RE_COVERAGE,
        readlen  = RE_READLEN,
        seed     = RE_SEED,
        workdir  = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "prep_{species}.log")
    resources:
        cpus_per_task  = 64,
        mem_mb_per_cpu = 2000,
        runtime        = 720
    envmodules:
        "SAMtools/1.21-GCC-13.3.0",
        "fastp/1.0.1-GCC-13.3.0"
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        S="{wildcards.species}"
        mkdir -p {params.workdir}
        cd {params.workdir}

        echo "[$(date)] === Genome size from .fai ==="
        # samtools faidx is idempotent; sum column 2 (sequence lengths).
        samtools faidx {input.genome}
        GENOME_SIZE=$(awk '{{sum+=$2}} END{{print sum}}' {input.genome}.fai)
        echo "Genome size: $GENOME_SIZE bp"

        echo "[$(date)] === QC (fastp; quality cutoff 10, >=95% above cutoff) ==="
        # Quality FILTER (not length-trim); uniform length is enforced later by the
        # hard-trim to {params.readlen}. Mirrors the RepeatExplorer preprocessing recipe.
        fastp -i {input.r1} -I {input.r2} \
          -o ${{S}}_R1.qc.fq.gz -O ${{S}}_R2.qc.fq.gz \
          --qualified_quality_phred 10 \
          --unqualified_percent_limit 5 \
          --n_base_limit 0 \
          --length_required {params.readlen} \
          --disable_adapter_trimming \
          --thread {resources.cpus_per_task} \
          --json ${{S}}_fastp.json --html ${{S}}_fastp.html

        SURVIVED=$(zcat ${{S}}_R1.qc.fq.gz | awk 'END{{print NR/4}}')

        # Pairs needed for target coverage: pairs = (coverage * genome) / (2 * readlen)
        N_PAIRS=$(awk -v c={params.coverage} -v g="$GENOME_SIZE" -v l={params.readlen} \
          'BEGIN{{printf "%d", (c*g)/(2*l)}}')

        echo "Coverage: {params.coverage}x | read length: {params.readlen} bp"
        echo "Pairs needed: $N_PAIRS | pairs surviving QC: $SURVIVED"

        # --- guards: never call seqtk with an empty/zero count ---
        [[ "$N_PAIRS" =~ ^[0-9]+$ ]] && [ "$N_PAIRS" -gt 0 ] || {{
          echo "ERROR: N_PAIRS invalid ([$N_PAIRS]). Check genome size/coverage/readlen." >&2
          exit 1
        }}
        [[ "$SURVIVED" =~ ^[0-9]+$ ]] && [ "$SURVIVED" -gt 0 ] || {{
          echo "ERROR: no reads survived QC (SURVIVED=[$SURVIVED]). Check fastp output." >&2
          exit 1
        }}
        if [ "$SURVIVED" -lt "$N_PAIRS" ]; then
          echo "WARNING: fewer pairs after QC ($SURVIVED) than needed for {params.coverage}x ($N_PAIRS)."
        fi

        echo "[$(date)] === Sample read PAIRS (shared seed -s {params.seed}) ==="
        # Identical seed + identical count on both mates recovers complete pairs.
        {params.seqtk} sample -s{params.seed} ${{S}}_R1.qc.fq.gz "$N_PAIRS" > ${{S}}_R1.sub.fq
        {params.seqtk} sample -s{params.seed} ${{S}}_R2.qc.fq.gz "$N_PAIRS" > ${{S}}_R2.sub.fq

        echo "[$(date)] === Hard-trim to {params.readlen} bp (enforce uniform length) ==="
        # trimfq -L drops reads shorter than L, guaranteeing uniform length.
        {params.seqtk} trimfq -L {params.readlen} ${{S}}_R1.sub.fq > ${{S}}_R1.trim.fq
        {params.seqtk} trimfq -L {params.readlen} ${{S}}_R2.sub.fq > ${{S}}_R2.trim.fq

        echo "[$(date)] === Convert to FASTA + rename with coded names ==="
        {params.seqtk} seq -A ${{S}}_R1.trim.fq \
          | awk 'BEGIN{{n=0}} /^>/{{n++; printf ">%07d/1\n", n; next}} {{print}}' > ${{S}}_R1.fa
        {params.seqtk} seq -A ${{S}}_R2.trim.fq \
          | awk 'BEGIN{{n=0}} /^>/{{n++; printf ">%07d/2\n", n; next}} {{print}}' > ${{S}}_R2.fa

        echo "[$(date)] === Interleave mates (drop broken pairs) ==="
        # Both files carry the same coded IDs in the same order, so mergepe keeps
        # pairing and alternates /1 and /2 as RepeatExplorer2 requires.
        {params.seqtk} mergepe ${{S}}_R1.fa ${{S}}_R2.fa > {output.interleaved}

        TOTAL=$(grep -c '^>' {output.interleaved})
        echo "Total reads (should be ~2 x $N_PAIRS): $TOTAL"
        echo "[$(date)] Done. Interleaved FASTA: {output.interleaved}"
        """


rule REPEATEXPLORER:
    """Run RepeatExplorer2 (TAREAN) on the interleaved reads via apptainer."""
    input:
        interleaved = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_reads_interleaved.fasta"),
    output:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer.done"),
    params:
        sif     = REPEATEXPLORER_SIF,
        taxon   = RE_TAXON,
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "repeatexplorer_{species}.log")
    resources:
        cpus_per_task  = 64,
        mem_mb_per_cpu = 2000,
        runtime        = 2880
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        export PYTHONHASHSEED=0
        S="{wildcards.species}"
        cd {params.workdir}

        # seqclust creates the output dir itself; remove any stale one first.
        rm -rf re_output
        echo "[$(date)] === RepeatExplorer2 ==="
        apptainer exec --bind "$PWD":/data/ {params.sif} \
          /repex_tarean/seqclust -p -c {resources.cpus_per_task} -A -r 100000000 \
          -tax {params.taxon} -m 0.001 \
          -v /data/re_output /data/${{S}}_reads_interleaved.fasta
        touch {output.done}
        echo "[$(date)] Done. Output in {params.workdir}/re_output"
        """
