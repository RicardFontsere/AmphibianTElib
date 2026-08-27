rule priority_all:
    """Convenience target for species that already have a TEtrimmer library:
    the linearised lib2.0 FASTA, the RepeatMasker run against it and its
    divergence landscape, and the priority table. Run `snakemake priority_all
    ...` to avoid pulling the rest of `rule all` (which would try to build
    TEtrimmer etc. for every species)."""
    input:
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority", "final_priority.table.tab"),
               species=SPECIES_WITH_TETRIMMER),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}_lib2.0.fa"),
               species=SPECIES_WITH_TETRIMMER),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.align.divsum"),
               species=SPECIES_WITH_TETRIMMER),
        expand(os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.ATELIB.divsum.html"),
               species=SPECIES_WITH_TETRIMMER)



rule PRIORITY_TABLE:
    """Post-process the TEtrimmer library into a per-family priority table:
    cd-hit-est dereplication, consensus length, genome copy-number, RepBase
    novelty screen (against a pre-built RepBase blast db) and Pfam domains,
    merged into final_priority.table.tab. Runs in place in TEtrimmer/priority/."""
    input:
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "TEtrimmer_consensus_merged.fasta"),
        genome  = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
    output:
        table = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority", "final_priority.table.tab"),
    params:
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "priority"),
        cdddir  = config["CDD_DIR"],
        repbase = config["REPBASE_DB"],
        script  = os.path.join(workflow.basedir, "scripts", "priority_table.sh"),
    log:
        os.path.join(LOG_DIR, "PriorityTable", "priority_{species}.log")
    resources:
        cpus_per_task  = 10,
        mem_mb_per_cpu = 8000,
        runtime        = 200
    envmodules:
        "CD-HIT/4.8.1-GCC-13.3.0",
        "BLAST+/2.16.0-gompi-2024a",
        "GCC/13.3.0"
    shell:
        """
        mkdir -p {params.workdir} $(dirname {log})
        cd {params.workdir}
        bash {params.script} \
            {input.library} \
            {input.genome} \
            {params.cdddir} \
            {params.repbase} \
            {resources.cpus_per_task} \
            &> {log}
        """


rule TETRIMMER_LIB2:
    """Linearise the TEtrimmer consensus library into the final `lib2.0` FASTA:
    every record is folded onto a single sequence line and each header is
    truncated at the first whitespace, so downstream tools see clean IDs.
    Output is named after the species (genus_epithet), e.g.
    Bufotes_viridis -> Bufotes_viridis_lib2.0.fa."""
    input:
        done    = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}.tetrimmer.done"),
        library = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "TEtrimmer_consensus_merged.fasta"),
    output:
        lib = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}_lib2.0.fa"),
    log:
        os.path.join(LOG_DIR, "PriorityTable", "lib2.0_{species}.log")
    resources:
        cpus_per_task  = 1,
        mem_mb_per_cpu = 4000,
        runtime        = 30
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.lib})
        awk '/^>/ {{sub(/[[:space:]].*/, ""); printf("\\n%s\\n", $0); next}} {{printf("%s", $0)}} END {{printf("\\n")}}' \
            < {input.library} \
          | awk 'NR > 1' > {output.lib} 2> {log}
        """


rule RMSK_LIB2:
    """Mask the genome with the freshly built lib2.0 library.

    Only the alignment (.align) and summary table (.tbl) are wanted downstream,
    so `-gff` / `-xsmall` are dropped from the command and the bulky by-products
    (.cat.gz, .masked) are removed once RepeatMasker finishes -- RepeatMasker has
    no switch to suppress them. `-gccalc` stays: the GC-corrected divergences it
    computes are what the Kimura landscape is built from."""
    input:
        genome = os.path.join(GENOMES_DIR_DONE, "{species}", "{species}_headers.fna"),
        lib    = os.path.join(GENOMES_DIR_DONE, "{species}", "TEtrimmer", "{species}_lib2.0.fa"),
    output:
        align = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.align"),
        tbl   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.tbl"),
        out   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.out"),
    params:
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK"),
    log:
        os.path.join(LOG_DIR, "RMSK", "rmsk_{species}.log")
    resources:
        cpus_per_task  = 64,
        mem_mb_per_cpu = 3200,
        runtime        = 1740
    envmodules:
        "RepeatMasker/4.2.3-foss-2025a",
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        S={wildcards.species}
        mkdir -p {params.workdir}
        cd {params.workdir}
        echo "[$S] RMSK_LIB2 $(date +%F_%T)"

        # rmblast uses 4 cores per parallel job, so -pa is cpus/4.
        RepeatMasker \
            -pa $(( {resources.cpus_per_task} / 4 )) \
            -a \
            -no_is \
            -gccalc \
            -s \
            -dir {params.workdir} \
            -lib {input.lib} \
            {input.genome}

        # Drop what we do not use; .align/.tbl/.out are the keepers.
        G=$(basename {input.genome})
        rm -f "$G".cat.gz "$G".cat "$G".masked "$G".out.gff
        echo "[$S] RMSK_LIB2 done $(date +%F_%T)"
        """


rule REPEAT_LANDSCAPE:
    """Kimura divergence landscape from the RMSK_LIB2 alignments:
    calcDivergenceFromAlign.pl (dfam-tetools container) summarises the .align
    into a .divsum, then createRepeatLandscapNoFamilies.pl renders it against the
    genome size pulled from the .tbl. The re-annotated alignment that
    calcDivergenceFromAlign.pl writes alongside the summary is temp() -- only the
    .divsum is kept (it feeds the priority table later)."""
    input:
        align = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.align"),
        tbl   = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}_headers.fna.tbl"),
    output:
        divsum         = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.align.divsum"),
        html           = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.ATELIB.divsum.html"),
        align_with_div = temp(os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK", "{species}.align_with_div")),
    params:
        sif     = DFAM_TETOOLS_SIF,
        script  = os.path.join(workflow.basedir, "scripts", "createRepeatLandscapNoFamilies.pl"),
        workdir = os.path.join(GENOMES_DIR_DONE, "{species}", "RMSK"),
    log:
        os.path.join(LOG_DIR, "RMSK", "landscape_{species}.log")
    resources:
        cpus_per_task  = 2,
        mem_mb_per_cpu = 4000,
        runtime        = 240
    envmodules:
        "BioPerl/1.7.8-GCCcore-11.3.0",
    shell:
        r"""
        exec &> {log}
        set -euo pipefail
        S={wildcards.species}
        test -s {params.script} || {{ echo "missing {params.script}"; exit 1; }}
        cd {params.workdir}
        echo "[$S] REPEAT_LANDSCAPE $(date +%F_%T)"

        apptainer exec --bind "$PWD":"$PWD" {params.sif} \
            calcDivergenceFromAlign.pl \
                -s {output.divsum} \
                -a {output.align_with_div} \
                {input.align}

        # Genome size for the landscape's y-axis normalisation.
        SIZE=$(grep -i 'total length' {input.tbl} | tr -s ' ' | cut -f 3 -d ' ')
        echo "[$S] genome size from .tbl: $SIZE bp"

        perl {params.script} -div {output.divsum} -g "$SIZE" > {output.html}
        echo "[$S] REPEAT_LANDSCAPE done $(date +%F_%T)"
        """
