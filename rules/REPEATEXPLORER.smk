def repeatexplorer_reads(wildcards):
    """Return the paired-end read pair (R1, R2) for a species.

    Species without a usable pair are filtered out of SPECIES_WITH_READS before
    the DAG is built, so this always finds a pair for the species it is called on.
    """
    r1, r2 = find_reads(wildcards.species)
    return {"r1": r1, "r2": r2}


rule SHORT_READ_PREP:
    """fastp QC, subsample to target coverage, hard-trim to uniform length,
    recode headers and interleave the mates for RepeatExplorer2."""
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
        cpus_per_task = 64, mem_mb_per_cpu = 2000, runtime = 720
    envmodules:
        "SAMtools/1.21-GCC-13.3.0", "fastp/1.0.1-GCC-13.3.0"
    shell:
        r"""
        exec &> {log}; set -euo pipefail
        S={wildcards.species}; cd {params.workdir}
        samtools faidx {input.genome}
        N=$(awk -v c={params.coverage} -v l={params.readlen} '{{g+=$2}} END{{printf "%d",(c*g)/(2*l)}}' {input.genome}.fai)
        fastp -i {input.r1} -I {input.r2} -o $S.qc_1.fq.gz -O $S.qc_2.fq.gz \
          --qualified_quality_phred 10 --unqualified_percent_limit 5 --n_base_limit 0 \
          --length_required {params.readlen} --disable_adapter_trimming --thread {resources.cpus_per_task} \
          -j $S.fastp.json -h $S.fastp.html
        {params.seqtk} sample -s{params.seed} $S.qc_1.fq.gz $N | {params.seqtk} trimfq -L {params.readlen} - \
          | {params.seqtk} seq -A - | awk '/^>/{{printf ">%07d/1\n",++n;next}}1' > $S.1.fa
        {params.seqtk} sample -s{params.seed} $S.qc_2.fq.gz $N | {params.seqtk} trimfq -L {params.readlen} - \
          | {params.seqtk} seq -A - | awk '/^>/{{printf ">%07d/2\n",++n;next}}1' > $S.2.fa
        {params.seqtk} mergepe $S.1.fa $S.2.fa > {output.interleaved}
        """


rule REPEATEXPLORER:
    """Run RepeatExplorer2 (TAREAN) on the interleaved reads."""
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
        cpus_per_task = 64, mem_mb_per_cpu = 2000, runtime = 2880
    shell:
        r"""
        exec &> {log}; set -euo pipefail
        export PYTHONHASHSEED=0; cd {params.workdir}; rm -rf re_output
        apptainer exec --bind "$PWD":/data/ {params.sif} \
          /repex_tarean/seqclust -p -c {resources.cpus_per_task} -A -r 100000000 \
          -tax {params.taxon} -m 0.001 -v /data/re_output /data/{wildcards.species}_reads_interleaved.fasta
        touch {output.done}
        """
