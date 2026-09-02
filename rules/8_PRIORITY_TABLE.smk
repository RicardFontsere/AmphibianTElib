# ===========================================================================
# Priority table
#
#   PRIORITY_TABLE -> per-family curation priority table from the TEtrimmer
#                     library + the RepeatMasker .divsum / .out
#
# NOTE on how this rule hooks onto TEtrimmer:
#   TEtrimmer_consensus_merged.fasta is deliberately NOT declared as an
#   input/output anywhere. TEtrimmer can exit non-zero after having written a
#   perfectly good library, and Snakemake deletes declared outputs of a failed
#   job -- we would lose the library. So TETRIMMER is tracked by the
#   {species}.tetrimmer.done sentinel it touches itself, downstream rules take
#   that sentinel as input (this is what makes the DAG edge), and the library
#   itself is passed as a params and checked by hand in the shell body.
#   REPEATMASKER does exactly the same.
# ===========================================================================


rule priority_all:
    """Convenience target: build ONLY the priority tables for species that already
    have a TEtrimmer library. Run `snakemake priority_all ...` to avoid pulling the
    rest of `rule all` (which would try to build TEtrimmer etc. for every species)."""
    input:
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "Priority", "final_priority.table.tab"),
               species=SPECIES_WITH_TETRIMMER)



rule PRIORITY_TABLE:
    """Post-process the TEtrimmer library into a per-family priority table:
    cd-hit-est dereplication, consensus length, genome occupancy and Kimura age
    from the RepeatMasker .divsum, defragmented/full-length copy counts from the
    RepeatMasker .out, RepBase novelty screen, CDD domains and a TRF satellite
    screen, merged into final_priority.table.tab. Runs in place in Priority/,
    at the same depth as RMDB/ and RMSK/."""
    input:
        # The TEtrimmer sentinel, not the library itself -- see the note at the
        # top of this file. This is the edge TETRIMMER -> PRIORITY_TABLE.
        done    = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done"),
        divsum  = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.align.divsum"),
        rmout   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.out"),
    output:
        table = os.path.join(GENOMES_DIR_DONE, "{species}", "Priority", "final_priority.table.tab"),
    params:
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "TEtrimmer_consensus_merged.fasta"),
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "Priority"),
        cdddir  = config["CDD_DIR"],
        repbase = config["REPBASE_DB"],
        script  = os.path.join(workflow.basedir, "scripts", "priority_table_withTRF.sh"),
    log:
        os.path.join(LOG_DIR, "PriorityTable", "priority_{species}.log")
    resources:
        cpus_per_task  = 10,
        mem_mb_per_cpu = 8000,
        runtime        = 600
    envmodules:
        "CD-HIT/4.8.1-GCC-13.3.0",
        "BLAST+/2.16.0-gompi-2024a",
        "GCC/13.3.0",
        "TRF/4.09.1-GCCcore-13.3.0"
    shell:
        r"""
        mkdir -p {params.workdir} $(dirname {log})
        exec &> {log}
        set -euo pipefail
        echo "[{wildcards.species}] PRIORITY_TABLE $(date +%F_%T)"

        # The library is a params, so Snakemake has not checked it exists.
        # TEtrimmer touches its .done sentinel even when it exits non-zero, so a
        # missing/empty library here means the TEtrimmer run really did fail.
        if [ ! -s "{params.library}" ]; then
            echo "[{wildcards.species}] ERROR: TEtrimmer library missing or empty: {params.library}" >&2
            echo "[{wildcards.species}] check the TEtrimmer log before re-running this rule." >&2
            exit 1
        fi

        cd {params.workdir}
        bash {params.script} \
            {params.library} \
            {input.divsum} \
            {input.rmout} \
            {params.cdddir} \
            {params.repbase} \
            {resources.cpus_per_task}

        echo "[{wildcards.species}] PRIORITY_TABLE done $(date +%F_%T) -> Priority/"
        """
