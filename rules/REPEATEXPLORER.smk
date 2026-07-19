def repeatexplorer_reads(wildcards):
    """Return the paired-end read pair (R1, R2) for a species.

    Species without a usable pair are filtered out of SPECIES_WITH_READS before
    the DAG is built, so this always finds a pair for the species it is called on.
    """
    r1, r2 = find_reads(wildcards.species)
    return {"r1": r1, "r2": r2}


rule SHORT_READ_PREP:
    """fastp QC, subsample to target coverage, hard-trim to uniform length,
    recode headers and interleave the mates for RepeatExplorer2.

    The QC'd (full-depth) fastqs are kept as outputs so the satMiner-style
    round 2 can re-subsample from them without re-running fastp."""
    input:
        unpack(repeatexplorer_reads),
        genome = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
    output:
        interleaved = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_reads_interleaved.fasta"),
        qc1 = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.qc_1.fq.gz"),
        qc2 = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.qc_2.fq.gz"),
    params:
        seqtk    = SEQTK,
        coverage = RE_COVERAGE,
        readlen  = RE_READLEN,
        seed     = RE_SEED,
        workdir  = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "prep_{species}.log")
    resources:
        cpus_per_task = 16, mem_mb_per_cpu = 1000, runtime = 250
    envmodules:
        "SAMtools/1.21-GCC-13.3.0",
        "fastp/1.0.1-GCC-13.3.0"
    shell:
        r"""
        exec &> {log}; set -euo pipefail
        S={wildcards.species}; cd {params.workdir}
        echo "[$S] === SHORT_READ_PREP $(date) ==="

        samtools faidx {input.genome}
        GENOME=$(awk '{{g+=$2}} END{{print g}}' {input.genome}.fai)
        N=$(awk -v c={params.coverage} -v l={params.readlen} -v g="$GENOME" 'BEGIN{{printf "%d",(c*g)/(2*l)}}')
        IN=$(zcat {input.r1} | awk 'END{{print NR/4}}')
        echo "[$S] genome length       : $GENOME bp"
        echo "[$S] input read pairs     : $IN"
        echo "[$S] pairs needed ({params.coverage}x): $N"

        fastp -i {input.r1} -I {input.r2} -o {output.qc1} -O {output.qc2} \
          --qualified_quality_phred 10 --unqualified_percent_limit 5 --n_base_limit 0 \
          --length_required {params.readlen} --thread {resources.cpus_per_task} \
          -j $S.fastp.json -h $S.fastp.html
        QC=$(zcat {output.qc1} | awk 'END{{print NR/4}}')
        echo "[$S] pairs surviving QC   : $QC  (dropped $((IN-QC)))"
        [ "$QC" -ge "$N" ] || echo "[$S] WARNING: only $QC pairs survived QC, fewer than the $N needed for {params.coverage}x"

        {params.seqtk} sample -s{params.seed} {output.qc1} $N | {params.seqtk} trimfq -L {params.readlen} - \
          | {params.seqtk} seq -A - | awk '/^>/{{printf ">%07d/1\n",++n;next}}1' > $S.1.fa
        {params.seqtk} sample -s{params.seed} {output.qc2} $N | {params.seqtk} trimfq -L {params.readlen} - \
          | {params.seqtk} seq -A - | awk '/^>/{{printf ">%07d/2\n",++n;next}}1' > $S.2.fa
        KEPT=$(grep -c '^>' $S.1.fa || true)
        echo "[$S] pairs after subsample+trim to {params.readlen} bp: $KEPT"

        {params.seqtk} mergepe $S.1.fa $S.2.fa > {output.interleaved}
        TOTAL=$(grep -c '^>' {output.interleaved} || true)
        echo "[$S] interleaved reads    : $TOTAL  (expected ~2 x $KEPT)"
        echo "[$S] === done $(date) ==="
        """


rule REPEATEXPLORER:
    """Run RepeatExplorer2 (TAREAN) on the round-1 interleaved reads."""
    input:
        interleaved = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_reads_interleaved.fasta"),
    output:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer.done"),
    params:
        sif = REPEATEXPLORER_SIF, taxon = RE_TAXON,
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "repeatexplorer_{species}.log")
    resources:
        cpus_per_task = 20, mem_mb_per_cpu = 12000, runtime = 3000
    shell:
        r"""
        exec &> {log}; set -euo pipefail
        export PYTHONHASHSEED=0; cd {params.workdir}; rm -rf re_output
        apptainer exec --bind "$PWD":/data/ {params.sif} \
          /repex_tarean/seqclust -p -c {resources.cpus_per_task} -A -r 204000000 \
          -tax {params.taxon} -m 0.001 -v /data/re_output /data/{wildcards.species}_reads_interleaved.fasta
        touch {output.done}
        """


# ---------------------------------------------------------------------------
# satMiner-style iterative mining (round 2): collect the repeats already found,
# deplete the reads that belong to them, and re-cluster the remainder so rarer
# satellites rise above RepeatExplorer's minimum cluster size.
# ---------------------------------------------------------------------------

rule COLLECT_KNOWN_REPEATS:
    """Build a 'known repeats' reference from the round-1 clusters.

    Ports satMiner's rexp_get_contigs_re2.py + extract_seq.py: per cluster, keep
    the most abundant contigs up to 50% of the cluster's cumulative coverage."""
    input:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer.done"),
    output:
        ref = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_known_repeats_r1.fasta"),
    params:
        # RE2 writes per-cluster contigs under .../clusters/dir_CL####/contigs.info.fasta
        clusters = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "re_output", "seqclust", "clustering", "clusters"),
        script   = os.path.join(workflow.basedir, "scripts", "collect_known_repeats.py"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "collect_known_{species}.log")
    resources:
        cpus_per_task = 1, mem_mb_per_cpu = 4000, runtime = 30
    shell:
        r"""
        exec &> {log}; set -euo pipefail
        python3 {params.script} {params.clusters} {output.ref}
        """


rule DEPLETE_READS:
    """Remove reads matching the known repeats, using bwa-mem2 + samtools.

    A larger (round-2) subsample is drawn from the QC reads and mapped to the
    known-repeat reference; only pairs whose BOTH mates stay unmapped are kept
    (samtools -f 12), then re-trimmed/recoded/interleaved exactly like round 1."""
    input:
        ref    = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_known_repeats_r1.fasta"),
        qc1    = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.qc_1.fq.gz"),
        qc2    = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.qc_2.fq.gz"),
        genome = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
    output:
        depleted = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_depleted_r2.fasta"),
    params:
        seqtk     = SEQTK,
        readlen   = RE_READLEN,
        seed      = RE_SEED,
        coverage2 = RE_COVERAGE_R2,
        workdir   = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "deplete_{species}.log")
    resources:
        cpus_per_task = 32, mem_mb_per_cpu = 4000, runtime = 720
    envmodules:
        "bwa-mem2/2.3-GCC-14.2.0", "SAMtools/1.22.1-GCC-14.2.0"
    shell:
        r"""
        exec &> {log}; set -euo pipefail
        S={wildcards.species}; cd {params.workdir}
        GENOME=$(awk '{{g+=$2}} END{{print g}}' {input.genome}.fai)
        N=$(awk -v c={params.coverage2} -v l={params.readlen} -v g="$GENOME" 'BEGIN{{printf "%d",(c*g)/(2*l)}}')
        echo "[$S] R2 target pairs ({params.coverage2}x): $N"

        # Enough-reads check: the round-2 subsample is drawn from the round-1 QC pool.
        QC=$(zcat {input.qc1} | awk 'END{{print NR/4}}')
        echo "[$S] QC pairs available   : $QC"
        [ "$QC" -ge "$N" ] || echo "[$S] WARNING: only $QC QC pairs available, fewer than the $N needed for {params.coverage2}x"

        # Larger subsample from the round-1 QC reads (shared seed keeps mates paired).
        {params.seqtk} sample -s{params.seed} {input.qc1} $N > r2_1.fq
        {params.seqtk} sample -s{params.seed} {input.qc2} $N > r2_2.fq

        # Map to the known repeats; -f 12 keeps records where the read AND its mate
        # are both unmapped, i.e. pairs that do NOT belong to already-found repeats.
        bwa-mem2 index {input.ref}
        bwa-mem2 mem -t {resources.cpus_per_task} {input.ref} r2_1.fq r2_2.fq \
          | samtools view -b -f 12 -F 256 - \
          | samtools collate -O -u - \
          | samtools fastq -1 clean_1.fq -2 clean_2.fq -0 /dev/null -s /dev/null -n
        echo "[$S] pairs kept after depletion: $(( $(wc -l < clean_1.fq) / 4 ))"

        # QC already enforced length >= readlen, so trimfq drops nothing here and the
        # two mate files stay in lockstep -> the coded IDs pair up as in round 1.
        {params.seqtk} trimfq -L {params.readlen} clean_1.fq | {params.seqtk} seq -A - \
          | awk '/^>/{{printf ">%07d/1\n",++n;next}}1' > d1.fa
        {params.seqtk} trimfq -L {params.readlen} clean_2.fq | {params.seqtk} seq -A - \
          | awk '/^>/{{printf ">%07d/2\n",++n;next}}1' > d2.fa
        {params.seqtk} mergepe d1.fa d2.fa > {output.depleted}
        echo "[$S] depleted interleaved reads: $(grep -c '^>' {output.depleted} || true)"
        """


rule REPEATEXPLORER_R2:
    """Run RepeatExplorer2 on the depleted reads (round 2) into re_output_r2."""
    input:
        depleted = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_depleted_r2.fasta"),
    output:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer_r2.done"),
    params:
        sif = REPEATEXPLORER_SIF, taxon = RE_TAXON,
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "repeatexplorer_r2_{species}.log")
    resources:
        cpus_per_task = 20, mem_mb_per_cpu = 12000, runtime = 3000
    shell:
        r"""
        exec &> {log}; set -euo pipefail
        export PYTHONHASHSEED=0; cd {params.workdir}; rm -rf re_output_r2
        apptainer exec --bind "$PWD":/data/ {params.sif} \
          /repex_tarean/seqclust -p -c {resources.cpus_per_task} -A -r 204000000 \
          -tax {params.taxon} -m 0.001 -v /data/re_output_r2 /data/{wildcards.species}_depleted_r2.fasta
        touch {output.done}
        """
