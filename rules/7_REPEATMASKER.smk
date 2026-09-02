# ===========================================================================
# RepeatMasker
#
#   REPEATMASKER    -> mask the genome with the TEtrimmer library
#   RMSK_DIVERGENCE -> Kimura divergence summary + repeat landscape
#
# Per species everything lives in GENOMES_DIR_DONE/{species}/RMSK/, i.e. at the
# same depth as RMDB/ and TEtrimmer/, not inside them.
# ===========================================================================


rule REPEATMASKER:
    input:
        genome  = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
        done = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done")
    output:
        align = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.align"),
        out   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.out"),
        tbl   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.tbl"),
    params:
        rmsk_dir = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK"),
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "TEtrimmer_consensus_merged.fasta"),
    log:
        os.path.join(LOG_DIR, "RepeatMasker", "repeatmasker_{species}.log")
    resources:
        cpus_per_task  = 32,          
        mem_mb_per_cpu = 2400,
        runtime        = 2200,        
    envmodules:
        "RepeatMasker/4.2.3-foss-2025a",
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        mkdir -p {params.rmsk_dir}
        cd {params.rmsk_dir}
        echo "[{wildcards.species}] REPEATMASKER $(date +%F_%T)"

        # -pa = cores / 4 (one rmblast job takes 4 cores)
        RepeatMasker \
          -pa 8 \
          -a \
          -no_is \
          -xsmall \
          -lib {params.library} \
          -dir {params.rmsk_dir} \
          {input.genome}

        echo "[{wildcards.species}] REPEATMASKER done $(date +%F_%T) -> RMSK/"
        """


rule RMSK_DIVERGENCE:
    """Kimura divergence summary (.divsum) + repeat landscape HTML from the
    RepeatMasker .align. """
    input:
        align = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.align"),
        tbl   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.tbl"),
    output:
        divsum = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.align.divsum"),
        html   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.divsum.html"),
    params:
        sif       = config["DFAM_TETOOLS_SIF"],
        landscape = os.path.join(workflow.basedir, "scripts", "createRepeatLandscapNoFamilies.pl"),
        rmsk_dir  = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK"),
    log:
        os.path.join(LOG_DIR, "RepeatMasker", "divergence_{species}.log")
    resources:
        cpus_per_task  = 2,
        mem_mb_per_cpu = 8000,
        runtime        = 600,
    envmodules:
        "BioPerl/1.7.8-GCCcore-14.2.0",
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        S={wildcards.species}
        cd {params.rmsk_dir}
        echo "[$S] RMSK_DIVERGENCE $(date +%F_%T)"

        apptainer exec --bind "$PWD":/data/ {params.sif} \
          calcDivergenceFromAlign.pl -s /data/$S.align.divsum /data/${{S}}_headers.fna.align

        # Genome size from the .tbl header line ("total length: <bp> bp (...)")
        SIZE=$(grep -i 'total length' {input.tbl} | tr -s ' ' | cut -f 3 -d ' ')
        echo "[$S] genome size $SIZE bp"
        perl {params.landscape} -div {output.divsum} -g $SIZE > {output.html}

        echo "[$S] RMSK_DIVERGENCE done $(date +%F_%T)"
        """
