def species_prefix(wildcards):
    """Build a short library prefix from a 'genus_epithet' species name.

    First 3 letters of the genus (lower-case) + first 3 letters of the
    epithet with its initial capitalised, e.g. 'Ranitomeya_variabilis' -> 'ranVar'.
    """
    genus, epithet = wildcards.species.split("_", 1)
    return genus[:3].lower() + epithet[:3].capitalize()


rule RENAME_RM_LIB:
    """Linearise the RepeatModeler library and rename its headers with a short species prefix."""
    input:
        rm_done  = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}.repeatmodeler.done"),
        families = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}-families.reclassified.fa")
    output:
        renamed = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}_rm1.0.fasta")
    params:
        prefix = species_prefix,
        script = os.path.join(workflow.basedir, "scripts", "renameRMDLconsensi.pl"),
        tmp    = os.path.join(GENOMES_DIR, "{species}", "RMDB", "{species}_rm1.0_temp.fasta")
    log:
        os.path.join(LOG_DIR, "RenameRM", "rename_{species}.log")
    envmodules:
        "Perl-bundle-CPAN/5.40.0-GCCcore-14.2.0"
    shell:
        """
        awk '/^>/ {{printf("\\n%s\\n",$0);next;}} {{printf("%s",$0);}} END {{printf("\\n");}}' \
            < {input.families} \
            | awk 'NR > 1' > {params.tmp} 2> {log}
        perl {params.script} {params.tmp} {params.prefix} {output.renamed} 2>> {log}
        rm -f {params.tmp}
        """
