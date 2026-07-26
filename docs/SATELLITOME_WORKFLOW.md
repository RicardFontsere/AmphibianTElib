# From raw reads to a satellitome library

How the RepeatExplorer2 branch of this pipeline turns unassembled Illumina reads into a
species-specific, named satDNA library that can be used for genome annotation and for
comparison between amphibian species.

The workflow follows the satMiner / satellitome protocol of Ruiz-Ruano et al.
(2016, *Sci Rep* 6:28333), with the classification thresholds as applied in the
*Chorthippus parallelus* (Genes 2023) and characid fish (Front Genet 2022) satellitomes.

```
reads ─┬─ SHORT_READ_PREP ── REPEATEXPLORER ─┬─ COLLECT_KNOWN_REPEATS ── DEPLETE_READS ── REPEATEXPLORER_R2 ─┐
       │                                     │                                                             │
       └─ FASTQC_RAW (reporting only)        └─────────────────────── SAT_COLLECT ◄────────────────────────┘
                                                                           │
   SAT_DIMERIZE ── SAT_HOMOLOGY ── SAT_GROUP ── SAT_QUANTIFY ── SAT_LIBRARY ──► {species}_satellitome_1.0.fasta
```

Rules 1-5 (`SHORT_READ_PREP` … `REPEATEXPLORER_R2`) live in `rules/6_REPEATEXPLORER.smk`,
rules 6-11 in `rules/7_SATELLITOME.smk`. Everything is per species; all outputs land under
`CompletedRuns/{species}/RepeatExplorer/` and `CompletedRuns/{species}/Satellitome/`.

---

## Step 1 — `SHORT_READ_PREP`: reads to a RepeatExplorer-shaped input

RepeatExplorer requires reads that are **uniform in length, high quality, adapter-free,
paired and interleaved**, at low coverage. Anything else biases the graph clustering.

* `fastp` — adapter detection, quality filter (Q10, ≤5% unqualified bases, **zero Ns** —
  RE2 discards reads with N anyway), and a hard trim/length filter to exactly `RE_READLEN`
  (101 bp). Reads shorter than that are dropped, not padded.
* `seqtk sample` — subsample to `RE_COVERAGE` (0.25×) of the genome, computed from the
  assembly length in the `.fai`. The same seed is used for both mates so pairs stay in
  lockstep.
* **New:** the sample is capped at `RE_MAX_READS` pairs (2 M). Coverage is the wrong knob
  for large genomes: `seqclust` does an all-to-all read comparison, so a 10 Gb amphibian
  genome at 0.25× (≈12 M pairs) will not finish. The log reports the effective coverage
  after capping — record it, it is the number that belongs in the methods.
* Read names are recoded to `>0000001/1` / `>0000001/2` and the mates interleaved
  (`seqtk mergepe`). RE2 needs the `/1` `/2` convention to treat the data as paired.

Kept: `{species}.qc_1.fq.gz`, `qc_2.fq.gz` (reused twice downstream),
`{species}_reads_interleaved.fasta`, fastp JSON/HTML.

## Step 2 — `REPEATEXPLORER`: round 1 clustering

`seqclust -p -A -c <cpus> -r <RAM kB> -m 0.01 -tax METAZOA3.0` in the RepeatExplorer2
container. Graph-based clustering of the read set; TAREAN runs on the clusters whose graph
topology looks like a tandem repeat and writes consensus monomers to
`TAREAN_consensus_rank_{1,2,3,4}.fasta` (1 = high confidence satellite, 2 = low confidence,
3 = putative LTR, 4 = other).

`-m 0.01` means a cluster must hold ≥0.01% of the reads to be reported — that is the
detection floor, and it is exactly why one round is not enough.

## Step 3 — `COLLECT_KNOWN_REPEATS`: what round 1 already explains

Python 3 port of satMiner's `rexp_get_contigs_re2.py`: for each cluster, contigs are ranked
by read depth and the top ones kept until they account for 50% of that cluster's coverage.
This is a compact stand-in for "the repeats we have already found".

## Step 4 — `DEPLETE_READS`: the satMiner trick

A fresh, larger sample (`RE_COVERAGE_R2` = 0.5×) is drawn from the same QC'd reads and
mapped to the known-repeat reference with `bwa-mem2`. Only pairs where **both mates are
unmapped** (`samtools view -f 12`) survive. The abundant satellites that dominated round 1
are removed from the read pool, so the rarer families that sat below the 0.01% floor now
rise above it. The depleted set is capped at `RE_MAX_READS` for the same reason as round 1.

## Step 5 — `REPEATEXPLORER_R2`: round 2 clustering

Identical `seqclust` call on the depleted reads, into `re_output_r2/`. Round-2 cluster sizes
are **not** comparable to round-1 sizes — they are fractions of a depleted read set. This is
the reason abundance has to be re-estimated later (step 10).

## Step 6 — `SAT_COLLECT`: gather the putative satDNAs

`scripts/collect_tarean_satellites.py` reads the TAREAN consensus files of both rounds,
keeps the ranks listed in `SAT_TAREAN_RANKS` (1,2 — ranks 3/4 are LTR elements and rDNA),
and writes one FASTA plus a metadata table (round, cluster, TAREAN rank/confidence, monomer
length, cluster size and genome proportion where `CLUSTER_TABLE.csv` provides them).
Sequence IDs are `R1_CL12_1_1` style so the round of origin is never lost.

## Step 7 — `SAT_DIMERIZE`: fix the circular permutation

`scripts/dimerator.py` (Python 3 port of Ruiz-Ruano's `dimerator.py`) tandem-repeats each
monomer. Two things are being fixed at once:

1. A monomer consensus is a **circular permutation** — a genomic array can start anywhere in
   it, so half of the real matches straddle the monomer boundary and are lost. A multimer
   restores the junction.
2. RepeatMasker cannot anchor an alignment in a sequence shorter than its seed. Short
   monomers are common (17-31 bp in this dataset) and simply vanish.

The refinement over the original script: instead of always doubling, each monomer is
repeated **until it reaches `SAT_MIN_MULTIMER_LEN` (200 nt), minimum two copies**. A 176 bp
monomer gets 2 copies; a 21 bp monomer gets 10. Doubling alone would leave the short ones
invisible — the same failure reported for a 21 nt satDNA in the 2025 Pyrgomorphidae/
Acrididae study, which only mapped once extended to pentamers.

## Step 8 — `SAT_HOMOLOGY`: all-against-all, one RepeatMasker run per sequence

`scripts/sat_homology.py` reimplements `rm_homology.py` in Python 3 (the original is Python 2
and needs RepeatMasker on `PATH`; that is what most of its GitHub issues are about).

RepeatMasker is run **once per sequence**, using that single sequence as the library and the
whole multimer set as the query. This is not a performance quirk — in a single pass with the
full library the entries compete for each region and only the best hit is reported, which
destroys exactly the pairwise similarities the classification depends on. Jobs are run in a
thread pool (`cpus/4`, since rmblast uses ~4 cores per job).

Output: a directional table of `identity = 100 − RepeatMasker divergence` plus the fraction
of the shorter monomer covered by the alignment.

## Step 9 — `SAT_GROUP`: variants, families, superfamilies

`scripts/sat_group.py` symmetrises the pairwise table (best of the two directions), drops
hits covering less than `SAT_MIN_COV` (50%) of the shorter monomer, and does single-linkage
clustering at three nested thresholds:

| identity | meaning | what happens |
|---|---|---|
| ≥ `SAT_VARIANT_ID` (95%) | same variant | collapsed — one consensus kept (the automated form of "manually selecting a consensus for sequences above 95% identity") |
| ≥ `SAT_FAMILY_ID` (80%) | variants of the same family | both kept in the library, same family number, letter suffix |
| ≥ `SAT_SUPERFAMILY_ID` (50%) | same superfamily | separate families, annotated as one superfamily |

The representative of a collapsed variant group is chosen by TAREAN rank, then cluster size,
then monomer length. Because the thresholds share one edge set, variants ⊂ families ⊂
superfamilies is guaranteed.

**The superfamily cutoff is the soft one in the literature.** The characid fish satellitome
uses >50%; the *Chorthippus* paper uses "any detectable homology below 80%" (set
`SAT_SUPERFAMILY_ID: 0` for that). The 95%/80% boundaries are stable across studies. State
whichever you adopt explicitly in the methods.

## Step 10 — `SAT_QUANTIFY`: abundance from reads, not from cluster sizes

Two independent random samples of `SAT_QUANT_READS` (5 M) reads are drawn from the QC'd
reads and RepeatMasked against the multimerised non-redundant library. Abundance per family
= masked bp / sampled bp, reported as mean ± sd over the two replicates, with mean
divergence per entry (the raw material for a repeat landscape).

This replaces cluster size for two reasons: round-2 clusters are computed on depleted reads,
and one family can be split over several clusters. A fixed read sample is what makes the
per-family percentages comparable **between species** — which is the whole point for a
comparative amphibian study.

Intervals are merged per read and per entry, so overlapping hits are not double counted; the
summary file also reports the total satDNA fraction of the genome per replicate.

## Step 11 — `SAT_LIBRARY`: naming and the deliverable

`scripts/name_satellitome.py` ranks families by summed abundance and applies the Ruiz-Ruano
convention: **species abbreviation + `Sat` + catalogue number + `-` + consensus monomer
length**, e.g. `RteSat01-176`. Families with several variants get a letter suffix
(`RteSat01-176A`, `RteSat01-180B`); single-variant families get none. The abbreviation is
derived as genus initial + first two letters of the epithet (`Rana_temporaria` → `Rte`),
overridable.

Deliverables in `CompletedRuns/{species}/Satellitome/`:

| file | use |
|---|---|
| `{species}_satellitome_1.0.fasta` | the library — one **monomer** per entry, `#Satellite` class, ordered by abundance |
| `{species}_satellitome_1.0.multimers.fasta` | the same library multimerised — **use this one with RepeatMasker** for genome annotation |
| `{species}_satellitome_1.0.tsv` | catalogue: name, family, variant, superfamily, monomer length, A+T%, abundance ± sd, divergence, TAREAN confidence, source cluster, collapsed IDs |
| `{species}_satellite_groups.tsv` | full variant/family/superfamily assignment, including collapsed sequences |
| `{species}_satellite_abundance.tsv`, `_quant_summary.tsv` | per-replicate abundances and total satDNA content |
| `{species}_satellite_pairwise.tsv` | the pairwise identity matrix behind the classification |

---

## Methods paragraph template

> Paired-end reads were quality-filtered and trimmed to a uniform length of 101 bp with
> fastp, and a random subset corresponding to ~X× genome coverage (N read pairs) was
> clustered with RepeatExplorer2 (`-p -A -tax METAZOA3.0 -m 0.01`). Reads matching the
> repeats identified in the first round (top contigs accounting for 50% of each cluster's
> coverage) were removed with bwa-mem2, and the remaining reads were re-clustered
> (satMiner-style iteration). TAREAN consensus sequences of high and low confidence
> satellite clusters from both rounds were tandem-multimerised to ≥200 nt (dimerator) and
> compared all-against-all with RepeatMasker, run sequence by sequence (rm_homology).
> Following Ruiz-Ruano et al. (2016), sequences sharing ≥95% identity were considered the
> same variant and collapsed into a single consensus, sequences with ≥80% identity were
> considered variants of the same family, and families with ≥50% identity were assigned to
> the same superfamily. Abundance was estimated by RepeatMasking two independent random
> samples of 5 million reads against the multimerised library and dividing masked bp by
> total sampled bp. Families were named as in Ruiz-Ruano et al. (2016): species
> abbreviation, "Sat", catalogue number in order of decreasing abundance, and consensus
> monomer length.

## Known limitations and suggested next steps

Not implemented here; listed in rough order of value for a multi-species amphibian survey.

1. **More than two rounds.** satMiner iterates until no new families appear, doubling the
   read input each round. Two rounds is a compromise; a third round on the depleted set is
   cheap because the input keeps shrinking. Worth adding as a loop with a stop criterion
   ("no new family above X% appeared"), driven by a `SAT_ROUNDS` config key.
2. **Cross-check against the TE library.** Some TAREAN clusters are tandemly arranged TE
   fragments or rDNA/organellar repeats. Masking the satellite library against the
   TEtrimmer/RepeatModeler consensi and flagging (not silently dropping) the overlaps would
   make the two branches of this pipeline consistent, and would stop the same sequence being
   annotated twice in the genome.
3. **Genomic validation on the assembly.** RepeatMasking the assembly with the multimerised
   library gives array lengths, chromosomal distribution and the assembled-vs-read abundance
   ratio — the usual sign of collapsed satellite arrays in an assembly, and useful in
   amphibians where the assemblies are far from complete.
4. **Tandem structure check.** TRF or a dotplot on the multimers filters out entries that
   are not really tandem, which TAREAN rank 2 occasionally lets through.
5. **Cross-species comparison.** The natural extension for the phylogenetic goal: run the
   same `sat_homology` + threshold logic **between** species libraries to define orthologous
   satDNA families across amphibians, then use presence/absence and divergence as characters.
   This is a new rule operating on the set of finished libraries rather than per species.
6. **Validation.** FISH remains the standard confirmation (e.g. the *Talpa aquitania*
   satellitome, Genes 2023); long reads (ONT/HiFi) are the practical alternative for
   confirming array organisation without cytogenetics.
