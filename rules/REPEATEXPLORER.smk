# ===========================================================================
# RepeatExplorer2 branch
#
#   round 1:  SHORT_READ_PREP -> REPEATEXPLORER
#   round 2:  COLLECT_KNOWN_REPEATS -> DEPLETE_READS -> REPEATEXPLORER_R2
#             (satMiner-style: remove reads of the repeats already found, then
#              re-cluster the remainder so rarer satellites rise above threshold)
#
# Per species everything lives in GENOMES_DIR_DONE/{species}/RepeatExplorer/.
# Kept outputs: qc_1/qc_2 fastqs (reused by round 2), the interleaved/depleted
# FASTAs (RE2 inputs), fastp reports, and re_output[_r2]/. Per-mate scratch
# files are removed at the end of each rule.
# ===========================================================================


def repeatexplorer_reads(wildcards):
    """R1/R2 read pair for a species (species without a pair are already
    filtered out of SPECIES_WITH_READS, so this always resolves)."""
    r1, r2 = find_reads(wildcards.species)
    return {"r1": r1, "r2": r2}


rule SHORT_READ_PREP:
    """Round 1 prep: fastp QC, subsample to target coverage, hard-trim to a
    uniform length, keep original read names, and interleave the mates."""
    input:
        unpack(repeatexplorer_reads),
        genome = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
    output:
        interleaved = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_reads_interleaved.fasta"),
        qc1         = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.qc_1.fq.gz"),
        qc2         = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.qc_2.fq.gz"),
        fastp_json  = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.fastp.json"),
        fastp_html  = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.fastp.html"),
    params:
        seqtk    = SEQTK,
        coverage = RE_COVERAGE,
        readlen  = RE_READLEN,
        seed     = RE_SEED,
        workdir  = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "prep_{species}.log")
    resources:
        cpus_per_task  = 16,
        mem_mb_per_cpu = 1000,
        runtime        = 250,
    envmodules:
        "SAMtools/1.21-GCC-13.3.0",
        "fastp/1.0.1-GCC-13.3.0",
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        S={wildcards.species}
        cd {params.workdir}
        echo "[$S] SHORT_READ_PREP $(date +%F_%T)"

        # Genome size (sum of .fai lengths) -> read pairs needed for the coverage.
        samtools faidx {input.genome}
        GENOME=$(awk '{{g+=$2}} END{{print g}}' {input.genome}.fai)
        N=$(awk -v c={params.coverage} -v l={params.readlen} -v g="$GENOME" 'BEGIN{{printf "%d",(c*g)/(2*l)}}')
        IN=$(zcat {input.r1} | awk 'END{{print NR/4}}')
        echo "[$S] genome $GENOME bp | target {params.coverage}x = $N pairs | input $IN pairs"

        # Quality FILTER only (uniform length is enforced below by trimfq -L).
        fastp -i {input.r1} -I {input.r2} -o {output.qc1} -O {output.qc2} \
          --qualified_quality_phred 10 --unqualified_percent_limit 5 --n_base_limit 0 \
          --length_required {params.readlen} --thread {resources.cpus_per_task} \
          -j {output.fastp_json} -h {output.fastp_html}
        QC=$(zcat {output.qc1} | awk 'END{{print NR/4}}')
        echo "[$S] QC: $QC pairs survived (dropped $((IN-QC)))"
        [ "$QC" -ge "$N" ] || echo "[$S] WARNING: QC pairs $QC < target $N"

        # Subsample -> trim to uniform length -> keep the original read name with a
        # /1 or /2 mate tag (so seqclust -k can preserve them) -> interleave mates.
        {params.seqtk} sample -s{params.seed} {output.qc1} $N | {params.seqtk} trimfq -L {params.readlen} - \
          | {params.seqtk} seq -A - | awk '/^>/{{print $1"/1"; next}}1' > $S.1.fa
        {params.seqtk} sample -s{params.seed} {output.qc2} $N | {params.seqtk} trimfq -L {params.readlen} - \
          | {params.seqtk} seq -A - | awk '/^>/{{print $1"/2"; next}}1' > $S.2.fa
        {params.seqtk} mergepe $S.1.fa $S.2.fa > {output.interleaved}
        echo "[$S] interleaved $(grep -c '^>' {output.interleaved}) reads ($(grep -c '^>' $S.1.fa) pairs)"

        rm -f $S.1.fa $S.2.fa   # drop per-mate scratch; keep qc_*/interleaved/reports
        """


rule REPEATEXPLORER:
    """Round 1 clustering: RepeatExplorer2 (TAREAN) on the interleaved reads."""
    input:
        interleaved = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_reads_interleaved.fasta"),
    output:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer.done"),
    params:
        sif     = REPEATEXPLORER_SIF,
        taxon   = RE_TAXON,
        mincl   = RE_MINCL,
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "repeatexplorer_{species}.log")
    resources:
        cpus_per_task  = 20,
        mem_mb_per_cpu = 12000,
        runtime        = 3000,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        export PYTHONHASHSEED=0
        cd {params.workdir}
        rm -rf re_output                       # seqclust must create the output dir itself
        RAM_KB=$(( {resources.cpus_per_task} * {resources.mem_mb_per_cpu} * 1000 ))
        # -p paired, -A auto-filter, -k keep read names, -r max RAM (kB), -m min cluster %
        apptainer exec --bind "$PWD":/data/ {params.sif} \
          /repex_tarean/seqclust -p -A -k -c {resources.cpus_per_task} -r $RAM_KB \
          -m {params.mincl} -tax {params.taxon} \
          -v /data/re_output /data/{wildcards.species}_reads_interleaved.fasta
        touch {output.done}
        """


rule COLLECT_KNOWN_REPEATS:
    """Round 2, step 1: build the 'known repeats' depletion reference from the
    round-1 clusters (top contigs to 50% of each cluster's coverage)."""
    input:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer.done"),
    output:
        ref = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_known_repeats_r1.fasta"),
    params:
        # per-cluster contigs live at .../clusters/dir_CL####/contigs.info.fasta
        clusters = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "re_output", "seqclust", "clustering", "clusters"),
        script   = os.path.join(workflow.basedir, "scripts", "collect_known_repeats.py"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "collect_known_{species}.log")
    resources:
        cpus_per_task  = 1,
        mem_mb_per_cpu = 4000,
        runtime        = 30,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        python3 {params.script} {params.clusters} {output.ref}
        """


rule DEPLETE_READS:
    """Round 2, step 2: subsample a larger read set from the round-1 QC reads,
    map it to the known repeats with bwa-mem2, and keep only the pairs that do
    NOT match (both mates unmapped) as the depleted, re-interleaved input."""
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
        cpus_per_task  = 32,
        mem_mb_per_cpu = 4000,
        runtime        = 720,
    envmodules:
        "bwa-mem2/2.3-GCC-14.2.0",
        "SAMtools/1.22.1-GCC-14.2.0",
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        S={wildcards.species}
        cd {params.workdir}
        echo "[$S] DEPLETE_READS $(date +%F_%T)"

        GENOME=$(awk '{{g+=$2}} END{{print g}}' {input.genome}.fai)
        N=$(awk -v c={params.coverage2} -v l={params.readlen} -v g="$GENOME" 'BEGIN{{printf "%d",(c*g)/(2*l)}}')
        QC=$(zcat {input.qc1} | awk 'END{{print NR/4}}')
        echo "[$S] genome $GENOME bp | R2 target {params.coverage2}x = $N pairs | QC pool $QC pairs"
        [ "$QC" -ge "$N" ] || echo "[$S] WARNING: QC pairs $QC < R2 target $N (using all available)"

        # Larger subsample from the round-1 QC reads (same seed keeps mates paired).
        {params.seqtk} sample -s{params.seed} {input.qc1} $N > r2_1.fq
        {params.seqtk} sample -s{params.seed} {input.qc2} $N > r2_2.fq
        SAMPLED=$(( $(wc -l < r2_1.fq) / 4 ))

        # Map to known repeats; -f 12 keeps records where the read AND its mate are
        # both unmapped, i.e. pairs that do NOT belong to already-found repeats.
        bwa-mem2 index {input.ref}
        bwa-mem2 mem -t {resources.cpus_per_task} {input.ref} r2_1.fq r2_2.fq \
          | samtools view -b -f 12 -F 256 - \
          | samtools collate -O -u - \
          | samtools fastq -1 clean_1.fq -2 clean_2.fq -0 /dev/null -s /dev/null -n
        KEPT=$(( $(wc -l < clean_1.fq) / 4 ))
        echo "[$S] depletion: $SAMPLED sampled -> $KEPT kept (removed $((SAMPLED-KEPT)) matching known repeats)"

        # QC already enforced length >= readlen, so trimfq drops nothing and the mates
        # stay in lockstep. Keep original names (+ /1,/2) as in round 1, then interleave.
        {params.seqtk} trimfq -L {params.readlen} clean_1.fq | {params.seqtk} seq -A - \
          | awk '/^>/{{print $1"/1"; next}}1' > d1.fa
        {params.seqtk} trimfq -L {params.readlen} clean_2.fq | {params.seqtk} seq -A - \
          | awk '/^>/{{print $1"/2"; next}}1' > d2.fa
        {params.seqtk} mergepe d1.fa d2.fa > {output.depleted}
        echo "[$S] depleted interleaved $(grep -c '^>' {output.depleted}) reads -> RepeatExplorer round 2"

        rm -f r2_1.fq r2_2.fq clean_1.fq clean_2.fq d1.fa d2.fa {input.ref}.*  # scratch + bwa index
        """


rule REPEATEXPLORER_R2:
    """Round 2, step 3: RepeatExplorer2 on the depleted reads into re_output_r2."""
    input:
        depleted = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}_depleted_r2.fasta"),
    output:
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer_r2.done"),
    params:
        sif     = REPEATEXPLORER_SIF,
        taxon   = RE_TAXON,
        mincl   = RE_MINCL,
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer"),
    log:
        os.path.join(LOG_DIR, "RepeatExplorer", "repeatexplorer_r2_{species}.log")
    resources:
        cpus_per_task  = 20,
        mem_mb_per_cpu = 12000,
        runtime        = 3000,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        export PYTHONHASHSEED=0
        cd {params.workdir}
        rm -rf re_output_r2
        RAM_KB=$(( {resources.cpus_per_task} * {resources.mem_mb_per_cpu} * 1000 ))
        apptainer exec --bind "$PWD":/data/ {params.sif} \
          /repex_tarean/seqclust -p -A -k -c {resources.cpus_per_task} -r $RAM_KB \
          -m {params.mincl} -tax {params.taxon} \
          -v /data/re_output_r2 /data/{wildcards.species}_depleted_r2.fasta
        touch {output.done}
        """
