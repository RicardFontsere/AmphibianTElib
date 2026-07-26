# ===========================================================================
# Satellitome branch: RepeatExplorer2/TAREAN output -> named satDNA library
#
#   SAT_COLLECT    gather TAREAN consensus monomers from both RE2 rounds
#   SAT_DIMERIZE   tandem-multimerise them to >= SAT_MIN_MULTIMER_LEN nt
#   SAT_HOMOLOGY   all-vs-all RepeatMasker, one run per sequence
#   SAT_GROUP      variants (>=95%) / families (>=80%) / superfamilies (>=50%)
#   SAT_QUANTIFY   abundance by RepeatMasking random read samples
#   SAT_LIBRARY    name by decreasing abundance -> {species}_satellitome_X.Y.fasta
#
# Everything lives in GENOMES_DIR_DONE/{species}/Satellitome/. The final
# deliverables are the named monomer library, its multimerised twin (the one
# to use with RepeatMasker for genome annotation) and the catalogue TSV.
#
# Classification criteria: Ruiz-Ruano et al. 2016 Sci Rep 6:28333; thresholds
# as applied in Ruiz-Ruano et al. 2023 (Genes 14:397) and Utsunomia et al.
# 2022 (Front Genet 13:891925).
# ===========================================================================


SATELLITOME_DIR = os.path.join(GENOMES_DIR_DONE, "{species}", "Satellitome")


rule SAT_COLLECT:
    """Pull the TAREAN consensus monomers of both RepeatExplorer2 rounds into a
    single putative-satDNA FASTA plus a metadata table."""
    input:
        r1_done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer.done"),
        r2_done = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.repeatexplorer_r2.done"),
    output:
        fasta = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.fasta"),
        table = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.tsv"),
    params:
        script = os.path.join(workflow.basedir, "scripts", "collect_tarean_satellites.py"),
        re1    = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "re_output"),
        re2    = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "re_output_r2"),
        ranks  = SAT_TAREAN_RANKS,
    log:
        os.path.join(LOG_DIR, "Satellitome", "Collect", "collect_{species}.log")
    resources:
        cpus_per_task  = 1,
        mem_mb_per_cpu = 4000,
        runtime        = 30,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        python3 {params.script} --out {output.fasta} --table {output.table} \
            --ranks {params.ranks} R1={params.re1} R2={params.re2}
        """


rule SAT_DIMERIZE:
    """Tandem-multimerise every monomer to at least SAT_MIN_MULTIMER_LEN nt.

    Monomers are circular permutations and short ones (17-31 bp is common) are
    invisible to RepeatMasker on their own; repeating them restores the junction
    and gives the aligner something to anchor on."""
    input:
        fasta = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.fasta"),
    output:
        dimers = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.multimers.fasta"),
        table  = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.multimers.tsv"),
    params:
        script  = os.path.join(workflow.basedir, "scripts", "dimerator.py"),
        min_len = SAT_MIN_MULTIMER_LEN,
    log:
        os.path.join(LOG_DIR, "Satellitome", "Dimerize", "dimerize_{species}.log")
    resources:
        cpus_per_task  = 1,
        mem_mb_per_cpu = 2000,
        runtime        = 15,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        python3 {params.script} --min-len {params.min_len} --table {output.table} \
            {input.fasta} {output.dimers}
        """


rule SAT_HOMOLOGY:
    """All-against-all homology: RepeatMasker run once per sequence (rm_homology
    logic), so related satDNAs do not lose their hit to a better-scoring library
    entry."""
    input:
        dimers = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.multimers.fasta"),
        table  = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.tsv"),
    output:
        pairwise = os.path.join(SATELLITOME_DIR, "{species}_satellite_pairwise.tsv"),
    params:
        script  = os.path.join(workflow.basedir, "scripts", "sat_homology.py"),
        cutoff  = SAT_RM_CUTOFF,
        workdir = os.path.join(SATELLITOME_DIR, "homology_work"),
    log:
        os.path.join(LOG_DIR, "Satellitome", "Homology", "homology_{species}.log")
    resources:
        cpus_per_task  = 16,
        mem_mb_per_cpu = 2000,
        runtime        = 600,
    envmodules:
        REPEATMASKER_MODULE,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        # rmblast uses ~4 cores per RepeatMasker job, so run cpus/4 jobs at once.
        JOBS=$(( {resources.cpus_per_task} / 4 )); [ "$JOBS" -ge 1 ] || JOBS=1
        rm -rf {params.workdir}
        python3 {params.script} --dimers {input.dimers} --lengths {input.table} \
            --out {output.pairwise} --threads $JOBS --cutoff {params.cutoff} \
            --workdir {params.workdir}
        rm -rf {params.workdir}
        """


rule SAT_GROUP:
    """Collapse redundancy and classify: >=95% identity = same variant (one
    consensus kept), >=80% = variants of one family, >=50% = one superfamily."""
    input:
        fasta    = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.fasta"),
        table    = os.path.join(SATELLITOME_DIR, "{species}_putative_satellites.tsv"),
        pairwise = os.path.join(SATELLITOME_DIR, "{species}_satellite_pairwise.tsv"),
    output:
        groups = os.path.join(SATELLITOME_DIR, "{species}_satellite_groups.tsv"),
        nr     = os.path.join(SATELLITOME_DIR, "{species}_satellites_nonredundant.fasta"),
        nr_dim = os.path.join(SATELLITOME_DIR, "{species}_satellites_nonredundant.multimers.fasta"),
    params:
        script    = os.path.join(workflow.basedir, "scripts", "sat_group.py"),
        dimerator = os.path.join(workflow.basedir, "scripts", "dimerator.py"),
        variant   = SAT_VARIANT_ID,
        family    = SAT_FAMILY_ID,
        superfam  = SAT_SUPERFAMILY_ID,
        min_cov   = SAT_MIN_COV,
        min_len   = SAT_MIN_MULTIMER_LEN,
    log:
        os.path.join(LOG_DIR, "Satellitome", "Group", "group_{species}.log")
    resources:
        cpus_per_task  = 1,
        mem_mb_per_cpu = 4000,
        runtime        = 30,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        python3 {params.script} --monomers {input.fasta} --table {input.table} \
            --pairwise {input.pairwise} --groups {output.groups} \
            --nonredundant {output.nr} \
            --variant-id {params.variant} --family-id {params.family} \
            --superfamily-id {params.superfam} --min-cov {params.min_cov}
        # Multimerised copy of the collapsed set: this is what the reads are
        # masked against in SAT_QUANTIFY.
        python3 {params.dimerator} --min-len {params.min_len} \
            --repeat-class Satellite {output.nr} {output.nr_dim}
        """


rule SAT_QUANTIFY:
    """Abundance re-estimation: RepeatMask two independent random read samples
    against the multimerised library and divide masked bp by sampled bp.

    Cluster sizes are not comparable across runs (round 2 is clustered on
    depleted reads); a fixed read sample is, which is what makes per-family
    percentages comparable between species."""
    input:
        nr_dim = os.path.join(SATELLITOME_DIR, "{species}_satellites_nonredundant.multimers.fasta"),
        groups = os.path.join(SATELLITOME_DIR, "{species}_satellite_groups.tsv"),
        qc1    = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.qc_1.fq.gz"),
        qc2    = os.path.join(GENOMES_DIR_DONE, "{species}", "RepeatExplorer", "{species}.qc_2.fq.gz"),
    output:
        abundance = os.path.join(SATELLITOME_DIR, "{species}_satellite_abundance.tsv"),
        summary   = os.path.join(SATELLITOME_DIR, "{species}_satellite_quant_summary.tsv"),
        out1      = os.path.join(SATELLITOME_DIR, "quant", "rep1.fa.out.gz"),
        out2      = os.path.join(SATELLITOME_DIR, "quant", "rep2.fa.out.gz"),
    params:
        script  = os.path.join(workflow.basedir, "scripts", "quantify_satellitome.py"),
        seqtk   = SEQTK,
        nreads  = SAT_QUANT_READS,
        seed    = RE_SEED,
        quantdir = os.path.join(SATELLITOME_DIR, "quant"),
    log:
        os.path.join(LOG_DIR, "Satellitome", "Quantify", "quantify_{species}.log")
    resources:
        cpus_per_task  = 20,
        mem_mb_per_cpu = 2000,
        runtime        = 1200,
    envmodules:
        REPEATMASKER_MODULE,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        S={wildcards.species}
        rm -rf {params.quantdir}; mkdir -p {params.quantdir}
        cd {params.quantdir}
        PA=$(( {resources.cpus_per_task} / 4 )); [ "$PA" -ge 1 ] || PA=1

        # Two independent random read samples (different seeds, one per mate file).
        {params.seqtk} sample -s$(( {params.seed} + 101 )) {input.qc1} {params.nreads} \
            | {params.seqtk} seq -A - > rep1.fa
        {params.seqtk} sample -s$(( {params.seed} + 202 )) {input.qc2} {params.nreads} \
            | {params.seqtk} seq -A - > rep2.fa

        printf 'replicate\tn_reads\ttotal_bp\n' > totals.tsv
        for R in rep1 rep2; do
            awk -v r=$R '/^>/{{n++; next}} {{b+=length($0)}} END{{printf "%s\t%d\t%d\n",r,n,b}}' \
                $R.fa >> totals.tsv
        done
        cat totals.tsv

        if [ -s {input.nr_dim} ]; then
            for R in rep1 rep2; do
                RepeatMasker -lib {input.nr_dim} -pa $PA -nolow -no_is -dir . $R.fa
            done
        else
            echo "[$S] WARNING: empty satDNA library, writing empty RepeatMasker output"
            for R in rep1 rep2; do
                printf '   SW  perc perc perc  query\n\n\n' > $R.fa.out
            done
        fi

        python3 {params.script} --totals totals.tsv --groups {input.groups} \
            --out {output.abundance} --summary {output.summary} \
            rep1=rep1.fa.out rep2=rep2.fa.out

        gzip -f rep1.fa.out rep2.fa.out
        rm -f rep1.fa rep2.fa *.masked *.cat.gz *.ori.out
        """


rule SAT_LIBRARY:
    """Final library: families numbered by decreasing abundance and named
    <Abb>Sat##-<monomer length>, variants suffixed A/B/..., plus the
    multimerised library used for genome annotation and the catalogue TSV."""
    input:
        nr        = os.path.join(SATELLITOME_DIR, "{species}_satellites_nonredundant.fasta"),
        groups    = os.path.join(SATELLITOME_DIR, "{species}_satellite_groups.tsv"),
        abundance = os.path.join(SATELLITOME_DIR, "{species}_satellite_abundance.tsv"),
    output:
        library   = os.path.join(SATELLITOME_DIR, "{species}_satellitome_" + SAT_LIB_VERSION + ".fasta"),
        multimers = os.path.join(SATELLITOME_DIR, "{species}_satellitome_" + SAT_LIB_VERSION + ".multimers.fasta"),
        catalogue = os.path.join(SATELLITOME_DIR, "{species}_satellitome_" + SAT_LIB_VERSION + ".tsv"),
    params:
        script    = os.path.join(workflow.basedir, "scripts", "name_satellitome.py"),
        dimerator = os.path.join(workflow.basedir, "scripts", "dimerator.py"),
        version   = SAT_LIB_VERSION,
        min_ab    = SAT_MIN_ABUNDANCE,
        min_len   = SAT_MIN_MULTIMER_LEN,
    log:
        os.path.join(LOG_DIR, "Satellitome", "Library", "library_{species}.log")
    resources:
        cpus_per_task  = 1,
        mem_mb_per_cpu = 2000,
        runtime        = 30,
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        python3 {params.script} --species {wildcards.species} \
            --nonredundant {input.nr} --groups {input.groups} \
            --abundance {input.abundance} --version {params.version} \
            --min-abundance {params.min_ab} \
            --out-fasta {output.library} --catalogue {output.catalogue}
        # Multimerised library for RepeatMasker-based genome annotation.
        python3 {params.dimerator} --min-len {params.min_len} \
            --repeat-class Satellite {output.library} {output.multimers}
        """
