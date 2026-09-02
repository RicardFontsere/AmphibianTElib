#!/bin/bash

set -e

if [ $# -ne 6 ]
then
    echo -e "\nusage: `basename $0` <TEtrimmer_lib.fa> <align.divsum> <tetrimmer_summary.txt> <cdd_dir> <repbase_blastdb> <threads>\n"
    echo -e "DESCRIPTION:\tRuns a pipeline that: 1) reduces sequence redundancy from the TEtrimmer library"
    echo -e "\t\tusing cd-hit-est; 2) extracts family names; 3) calculates consensus length; 4) pulls"
    echo -e "\t\twellCharLen and Kimura divergence from a RepeatMasker align.divsum; 4b) pulls"
    echo -e "\t\tBLAST copy counts from the TEtrimmer summary.txt; 5) screens each"
    echo -e "\t\tfamily against RepBase for novelty; 6) searches conserved domains against CDD;"
    echo -e "\t\t6b) screens each consensus with TRF and flags families whose length is >=80%"
    echo -e "\t\tcovered by a tandem array as satellite candidates; 7) merges into"
    echo -e "\t\tfinal_priority.table.tab.\n"
    echo -e "REQUIRES:\tcd-hit, trf, blast (incl. rpstblastn) and rpsbproc + a CDD rpsblast database.\n"
    echo -e "INPUT:\t\t<TEtrimmer_lib.fa>\tTEtrimmer consensus library (e.g. TEtrimmer_consensus_merged.fasta)"
    echo -e "\t\t<align.divsum>\t\tRepeatMasker .divsum file (Kimura table, class/repeat/absLen/wellCharLen/Kimura%)"
    echo -e "\t\t<tetrimmer_summary.txt>\tTEtrimmer summary.txt (CSV) from the run that built the library"
    echo -e "\t\t<cdd_dir>\t\tCDD base dir: holds the rpsbproc binary and utils/db/Cdd + utils/data"
    echo -e "\t\t<repbase_blastdb>\tpath to a PRE-BUILT RepBase blast nucl database (db prefix, already makeblastdb'd)"
    echo -e "\t\t<threads>\t\tnumber of CPU threads"
    echo -e "OUTPUT:\t\tfinal_priority.table.tab, with columns:"
    echo -e "\t\tfamily | length | input_full_blast_n | output_blast_n | output_full_blast_n | KimuraDiv |"
    echo -e "\t\trepbase_hit | repbase_pident | repbase_qcov | repbase_bitscore |"
    echo -e "\t\tn_domains | domains | trf_call | trf_copies | trf_monomer_len"
    echo -e "\t\tsat_monomers.nr.fa, dereplicated monomers of the flagged families\n"
    exit
fi

lib=$1
divsum=$2
summary=$3
cdddir=$4
repbase=$5
threads=$6


# ---------------------------------------------------------------------------
# P1  reduce redundancy
# ---------------------------------------------------------------------------

if [ ! -f cdhit.fa.clstr ]; then
    echo ">>> [P1] running cd-hit-est"
    cd-hit-est -i $lib -o cdhit.fa -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T $threads -M 0
fi

echo ">>> [P1] DONE - redundancy reduced"


# ---------------------------------------------------------------------------
# P2  family names  ->  col1.txt
# ---------------------------------------------------------------------------

awk '{ if ($1~/^>/) {print $1}}' cdhit.fa | sed 's/>//1' | sort > col1.txt

echo ">>> [P2] DONE - `wc -l < col1.txt` families after dereplication"


# ---------------------------------------------------------------------------
# P3  consensus length  ->  col2.txt
# ---------------------------------------------------------------------------

awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(substr(sequence_id, 2)" "sequence_length)}' cdhit.fa | grep "\S" > cdhit.fa.len
sort cdhit.fa.len | awk '{print $NF}' > col2.txt

echo ">>> [P3] DONE - consensus lengths calculated"


# ---------------------------------------------------------------------------
# P4  genome occupancy + age from align.divsum  ->  col_wcl.txt, col_kim.txt
# ---------------------------------------------------------------------------
# The divsum Kimura table is  Class <TAB> Repeat <TAB> absLen <TAB> wellCharLen <TAB> Kimura%.
# Keep only data rows (5 fields, numeric Kimura). If the same repeat name appears under
# more than one class, wellCharLen is summed and Kimura is averaged weighted by wellCharLen.
# Keys are the repeat name; library names are keyed on the part BEFORE '#' to match.

awk -F'\t' '
  NF>=5 && $5 ~ /^[0-9]+(\.[0-9]+)?$/ {
    w[$2]  += $4;
    kw[$2] += $5*$4;
  }
  END { for (r in w) printf "%s\t%d\t%.2f\n", r, w[r], (w[r]>0 ? kw[r]/w[r] : 0) }
' "$divsum" | sort > divsum.tab

awk -F'\t' '
  NR==FNR { wcl[$1]=$2; kim[$1]=$3; next }
  { k=$0; sub(/#.*/,"",k);
    print (k in wcl ? wcl[k] : "NA") > "col_wcl.txt";
    print (k in kim ? kim[k] : "NA") > "col_kim.txt";
  }
' divsum.tab col1.txt

echo ">>> [P4] DONE - `awk '$1!="NA"' col_wcl.txt | wc -l` families found in $divsum"


# ---------------------------------------------------------------------------
# P4b  BLAST copy counts from the TEtrimmer summary  ->  col_ifbn/col_obn/col_ofbn
# ---------------------------------------------------------------------------
# TEtrimmer already BLASTs every consensus back against the genome and records how
# many hits it recovered and how many of those were full length, for the input
# (RepeatModeler) consensus and for its own curated output. That answers the same
# question the RepeatMasker .out was being parsed for, measured by the tool that
# actually built the library, so the .out parsing is gone entirely -- along with its
# fragment defragmentation, its reconstructed consensus length, and the
# FULLLEN_MINCOV / COUNT_MINCOV thresholds those needed.
#
# summary.txt is CSV with a header row. Fields used:
#   2 output_name  4 input_full_blast_n  5 output_blast_n  6 output_full_blast_n
#
# Keyed on output_name -- the name TEtrimmer wrote into the library -- matched against
# the part of the library name BEFORE '#', the same keying P4/P6 use. One output_name
# can appear on several rows when TEtrimmer collapses inputs, and the output_* counts
# are properties of the output that repeat identically across those rows, so max()
# dereplicates them; summing would multiply-count a single consensus.
#
# TEtrimmer writes literal NaN for values it did not compute; those become NA so the
# column stays alignable with the rest of the table.

if [ ! -f "$summary" ]; then
    echo "ERROR: TEtrimmer summary not found: $summary" >&2
    exit 1
fi

awk -F',' '
  NR==1 { next }
  { sub(/\r$/, "") }
  $2 == "" || $2 == "NaN" { next }
  {
    k = $2; seen[k] = 1
    if ($4 != "NaN" && $4 != "" && (!(k in a) || $4+0 > a[k])) a[k] = $4+0
    if ($5 != "NaN" && $5 != "" && (!(k in b) || $5+0 > b[k])) b[k] = $5+0
    if ($6 != "NaN" && $6 != "" && (!(k in c) || $6+0 > c[k])) c[k] = $6+0
  }
  END {
    for (k in seen)
      printf "%s\t%s\t%s\t%s\n", k, (k in a ? a[k] : "NA"), \
                                        (k in b ? b[k] : "NA"), \
                                        (k in c ? c[k] : "NA")
  }
' "$summary" | sort > tetrimmer.tab

awk -F'\t' '
  NR==FNR { ifb[$1]=$2; obn[$1]=$3; ofbn[$1]=$4; next }
  { k=$0; sub(/#.*/,"",k);
    print (k in ifb ? ifb[k]  : "NA") > "col_ifbn.txt";
    print (k in ifb ? obn[k]  : "NA") > "col_obn.txt";
    print (k in ifb ? ofbn[k] : "NA") > "col_ofbn.txt";
  }
' tetrimmer.tab col1.txt

echo ">>> [P4b] DONE - `awk '$1!="NA"' col_obn.txt | wc -l` of `wc -l < col1.txt` families found in $summary"


# ---------------------------------------------------------------------------
# P5  novelty screen against RepBase  ->  col5.txt   (db is pre-built)
# ---------------------------------------------------------------------------

# Minimum HSP length (bp) for a RepBase hit to count. A short high-identity HSP
# (e.g. 95% over 40 bp) is low-complexity / nested / motif noise, not family
# evidence, so below this the family stays no-hit and goes to the novelty pool.
REPBASE_MINLEN=100

if [ ! -f repbase.blast ] || [ "$(awk 'NR==1{print NF; exit}' repbase.blast)" != "7" ]; then
    echo ">>> [P5] blasting library against pre-built RepBase db (dc-megablast, dust on)"
    blastn -task dc-megablast -query cdhit.fa -db $repbase -num_threads $threads -dust yes \
        -outfmt "6 qseqid sseqid pident length qcovs slen bitscore" -evalue 1e-5 > repbase.blast
fi

# Drop HSPs < MINLEN bp, rank remaining by BITSCORE, keep the best per family;
# record hit, pident, qcov and bitscore.
awk -v m="$REPBASE_MINLEN" '$4>=m' repbase.blast \
  | sort -k1,1 -k7,7gr \
  | awk '!s[$1]++ {printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$5,$7}' \
  | sort > repbase.best
join -a1 -e "0" -o 0,2.2,2.3,2.4,2.5 col1.txt repbase.best | tr ' ' '\t' | cut -f2-5 > col5.txt

echo ">>> [P5] DONE - `awk '$1!="0"' col5.txt | wc -l` families have a RepBase hit (>= ${REPBASE_MINLEN} bp)"


# ---------------------------------------------------------------------------
# P6  Conserved-domain search via CDD  ->  col4.txt (n_domains), col6.txt (domains)
# ---------------------------------------------------------------------------
# rpstblastn 6-frame-translates the nucleotide consensus and searches CDD PSSMs
# (Pfam + SMART + COG + PRK + TIGRFAM + curated CD) -- no ORF step, so degraded /
# frameshifted domains still hit. rpsbproc renders the readable hit table; we map
# its Query_N ids back to the real names, key on the family id BEFORE '#', and emit
# one line per col1 family (in col1 order) so col4/col6 stay aligned for the paste.

if [ ! -f cdd.asn ]; then
    echo ">>> [P6] rpstblastn vs CDD (6-frame)"
    rpstblastn -query cdhit.fa -db "$cdddir/utils/db/Cdd" -num_threads $threads -evalue 0.01 -outfmt 11 -out cdd.asn
fi

if [ ! -f cdd.proc ]; then
    echo ">>> [P6] rpsbproc post-processing"
    "$cdddir/rpsbproc" -i cdd.asn -o cdd.proc -d "$cdddir/utils/data" -m rep -e 0.01
fi

# Resolve Query_N -> real name, DROP Non-specific hits, emit  fam<TAB>from<TAB>domain;
# sort by family then genomic position so domains list N->C -- their order along the
# consensus is informative for classification (e.g. PR-RT-RH-IN in retroelements).
awk -F'\t' '
  $1=="QUERY"      { qname[$2]=$5; next }
  $1=="DOMAINS"    { indom=1; next }
  $1=="ENDDOMAINS" { indom=0; next }
  indom && $3!="Non-specific" {
    q=$2; sub(/\[.*/,"",q); fam=qname[q]; sub(/#.*/,"",fam);
    print fam"\t"$5"\t"$10
  }
' cdd.proc | sort -t$'\t' -k1,1 -k2,2n > cdd.domlist

# Aggregate per family in that positional order, one line per col1 family (aligned).
awk -F'\t' '
  NR==FNR { cnt[$1]++; nm[$1]=(nm[$1]==""?$3:nm[$1]","$3); next }
  { k=$0; sub(/#.*/,"",k);
    print (k in cnt ? cnt[k] : 0)      > "col4.txt";
    print (k in nm  ? nm[k]  : "none") > "col6.txt";
  }
' cdd.domlist col1.txt

echo ">>> [P6] DONE - CDD domain search finished"


# ---------------------------------------------------------------------------
# P7  tandem-repeat screen with TRF  ->  col_trf*.txt + sat_monomers.nr.fa
# ---------------------------------------------------------------------------
# maxperiod=2000 is the critical parameter: TRF's default of 500 silently misses
# the kb-scale monomers that are common in amphibians.
#
# -ngs stream fields per array:
#   1 start  2 end  3 period  4 copies  5 consensusSize  6 %match  7 %indel
#   8 score  9-12 A C G T  13 entropy  14 monomer  15 array_seq
#
# Coverage is the UNION of the reported start-end intervals, with the period left
# out of the key. TRF describes one array several times -- at the monomer period, at
# the dimer (HOR), and at drifted periods caused by indels -- so keying on the period
# splits a single satellite into fragments, and summing copies x monomerLen counts
# the same bases two or three times. The union counts every base once and can never
# exceed the consensus length, so trf_frac stays a real fraction.
#
# TRF does NOT emit arrays in ascending start order. It reports by increasing period,
# so the enclosing HOR description routinely arrives AFTER the narrower sub-arrays it
# contains -- on scaCou3_397_rnd3#rRNA the 9214-10822 array at period 67 is emitted
# last, after three 24-mers that sit inside it. A single pass that only tracks the
# running end then credits that HOR with its 41 bp tail instead of its 1609 bp span,
# giving a union of 1186 bp where the true answer is 2307. The error falls hardest on
# higher-order satellite structure, i.e. exactly what this screen is looking for, so
# the intervals are flattened and sorted before being merged:
#   start > e + 1        -> disjoint, bank the open run and start a new one
#   start <= e + 1       -> contiguous or overlapping, extend the run to max(e, end)
#
# TRF_MINPER=10 drops poly-A tails and microsatellites from the union; without it a
# real TE with a long homopolymer tail can be pushed over the threshold. Satellite
# monomers are well above 10 bp, so nothing genuine is lost.
#
# The exported monomer comes from the dominant array (most copies), tie-broken to the
# SMALLER period so the monomer is reported, not its dimer.
#
# A family is called SAT_CANDIDATE only when >= TRF_MINFRAC of its consensus length
# is covered; everything else is NA. CDD domains are deliberately NOT used to veto
# the call -- flagged families go to manual inspection, where a real TE carrying an
# internal tandem block (Helitron, hAT-Ac subterminal repeats, some TcMar) is
# separated from a satellite that TEtrimmer swallowed.

TRF_MINFRAC=0.8
TRF_MINPER=10

if [ ! -f trf.dat ]; then
    echo ">>> [P6b] running TRF (maxperiod 2000)"
    trf cdhit.fa 2 7 7 80 10 50 2000 -h -d -ngs > trf.dat || true
fi

paste col1.txt col2.txt > trf.keys.tsv

# Flatten the -ngs stream to  family / start / end / period / copies / consSize / monomer
# and sort by family then start, so the merge below really does see ascending starts.
awk -v OFS='\t' '
  /^@/   { id = substr($1,2); next }
  NF>=14 { print id, $1, $2, $3, $4, $5, $14 }
' trf.dat | sort -t$'\t' -k1,1 -k2,2n > trf.arrays.tsv

awk -F'\t' -v minfrac="$TRF_MINFRAC" -v minper="$TRF_MINPER" '
  NR==FNR { L[$1]=$2; next }
  $4+0 < minper { next }
  {
    if      ($1 != cur) { if (cur != "") cov[cur] += e-s+1; cur=$1; s=$2; e=$3 }
    else if ($2 > e+1)  { cov[cur] += e-s+1; s=$2; e=$3 }
    else if ($3 > e)    { e=$3 }
    # Dominant array: most copies, tie-broken to the SMALLER period so the monomer is
    # reported and not its dimer. A max over all arrays, so the sort does not affect it.
    if ($5+0 > n[$1] || ($5+0 == n[$1] && $4+0 < p[$1])) {
      n[$1]=$5; p[$1]=$4; ml[$1]=$6; mono[$1]=$7
    }
  }
  END {
    if (cur != "") cov[cur] += e-s+1          # bank the last open run
    # One line per family TRF found any array in. trf_frac is reported for all of
    # them so the near-misses stay visible; the copies / monomer_len / monomer
    # columns are filled only for families that clear the threshold.
    for (i in cov) {
      f = (L[i]>0 ? cov[i]/L[i] : 0)
      if (f >= minfrac)
        printf "%s\tSAT_CANDIDATE\t%.1f\t%s\t%.2f\t%s\n", i, n[i]+0, ml[i], f, mono[i]
      else
        printf "%s\tNA\tNA\tNA\t%.2f\tNA\n", i, f
    }
  }
' trf.keys.tsv trf.arrays.tsv | sort > trf.tab

awk -F'\t' '
  NR==FNR { cl[$1]=$2; c[$1]=$3; ml[$1]=$4; fr[$1]=$5; next }
  { print ($0 in cl ? cl[$0] : "NA") > "col_trfcall.txt";
    print ($0 in cl ? c[$0]  : "NA") > "col_trfcop.txt";
    print ($0 in cl ? ml[$0] : "NA") > "col_trfmonlen.txt";
    print ($0 in cl ? fr[$0] : "NA") > "col_trffrac.txt";
  }
' trf.tab col1.txt

# Monomers of the flagged families. BLAST-extend-extend pipelines enter the same
# satellite at several phases of the monomer, so clustering is not optional.
awk -F'\t' '$2=="SAT_CANDIDATE" {print ">"$1"_mono_L"$4"\n"$6}' trf.tab > sat_monomers.fa
if [ -s sat_monomers.fa ]; then
    cd-hit-est -i sat_monomers.fa -o sat_monomers.nr.fa -d 0 -c 0.8 -aS 0.8 -G 0 -T $threads -M 0 > /dev/null
fi

echo ">>> [P6b] DONE - `grep -c SAT_CANDIDATE col_trfcall.txt || true` families >= ${TRF_MINFRAC} covered by a tandem array"


# ---------------------------------------------------------------------------
# P8  merge
# ---------------------------------------------------------------------------

# paste silently produces a ragged table if any column file is short, so check first.
ncol=`wc -l < col1.txt`
for f in col2.txt col_ifbn.txt col_obn.txt col_ofbn.txt col_kim.txt col5.txt col4.txt col6.txt col_trfcall.txt col_trfcop.txt col_trfmonlen.txt; do
    n=`wc -l < $f`
    if [ "$n" -ne "$ncol" ]; then
        echo "ERROR: $f has $n lines, expected $ncol (same as col1.txt) - aborting merge" >&2
        exit 1
    fi
done

echo -e "family\tlength\tinput_full_blast_n\toutput_blast_n\toutput_full_blast_n\tKimuraDiv\trepbase_hit\trepbase_pident\trepbase_qcov\trepbase_bitscore\tn_domains\tdomains\ttrf_call\ttrf_copies\ttrf_monomer_len" > final_priority.table.tab
paste -d "\t" col1.txt col2.txt col_ifbn.txt col_obn.txt col_ofbn.txt col_kim.txt col5.txt col4.txt col6.txt col_trfcall.txt col_trfcop.txt col_trfmonlen.txt >> final_priority.table.tab

echo ">>> [P7] DONE - final_priority.table.tab generated"
echo
echo "To list candidate families with no close RepBase match, ranked by full-length BLAST hits:"
echo "  awk -F'\t' 'NR>1 && (\$7==0 || \$9<50)' final_priority.table.tab | sort -k5,5nr | less -S"
echo
echo "To list intact, young novel families (the best curation targets):"
echo "  awk -F'\t' 'NR>1 && (\$7==0 || \$9<50) && \$5>=5 && \$6!=\"NA\" && \$6<20' final_priority.table.tab | sort -k4,4nr | less -S"
echo
echo "To list satellite candidates for manual inspection, ranked by copy number:"
echo "  awk -F'\t' 'NR>1 && \$13==\"SAT_CANDIDATE\"' final_priority.table.tab | sort -k4,4nr | less -S"
echo
echo "Covered fraction is not in the table; near-misses are in trf.tab (field 5):"
echo "  awk -F'\t' '\$2==\"NA\" && \$5>=0.5' trf.tab | sort -k5,5gr | less -S"