#!/bin/bash

set -e

if [ $# -ne 6 ]
then
    echo -e "\nusage: `basename $0` <TEtrimmer_lib.fa> <genome.fa> <cdd_dir> <repbase_blastdb> <threads> <satellite_lib.fa>\n"
    echo -e "DESCRIPTION:\tRuns a pipeline that: 1) reduces sequence redundancy from the TEtrimmer library"
    echo -e "\t\tusing cd-hit-est; 2) extracts family names; 3) calculates consensus length; 4) makes a rough"
    echo -e "\t\testimate of genome copy number; 5) screens each family against RepBase for novelty;"
    echo -e "\t\t6) searches conserved domains against CDD; 7) merges into final_priority.table.tab.\n"
    echo -e "REQUIRES:\tcd-hit, blast (incl. rpstblastn) and rpsbproc + a CDD rpsblast database.\n"
    echo -e "INPUT:\t\t<TEtrimmer_lib.fa>\tTEtrimmer consensus library (e.g. TEtrimmer_consensus_merged.fasta)"
    echo -e "\t\t<genome.fa>\t\tgenome used to predict the library"
    echo -e "\t\t<cdd_dir>\t\tCDD base dir: holds the rpsbproc binary and utils/db/Cdd + utils/data"
    echo -e "\t\t<repbase_blastdb>\tpath to a PRE-BUILT RepBase blast nucl database (db prefix, already makeblastdb'd)"
    echo -e "\t\t<threads>\t\tnumber of CPU threads"
    echo -e "\t\t<satellite_lib.fa>\tsatellite consensus library FASTA (screened to flag satellite-only TEs)"
    echo -e "OUTPUT:\t\tfinal_priority.table.tab, with columns:"
    echo -e "\t\tfamily | length | copies | n_domains | repbase_hit | repbase_pident | repbase_qcov | repbase_alnlen | repbase_scov | repbase_bitscore | domains | sat_hit | sat_qcov\n"
    exit
fi

rmout=$1
genome=$2
cdddir=$3
repbase=$4
threads=$5
satlib=$6


# ---------------------------------------------------------------------------
# P1  reduce redundancy
# ---------------------------------------------------------------------------

if [ ! -f cdhit.fa.clstr ]; then
    echo ">>> [P1] running cd-hit-est"
    cd-hit-est -i $rmout -o cdhit.fa -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T $threads -M 0
fi

echo ">>> [P1] DONE - redundancy reduced"


# ---------------------------------------------------------------------------
# P2  family names  ->  col1.txt
# ---------------------------------------------------------------------------

awk '{ if ($1~/^>/) {print $1}}' cdhit.fa | sed 's/>//1' | sort > col1.txt

echo ">>> [P2] DONE - `wc -l < col1.txt` families after dereplication"


# ---------------------------------------------------------------------------
# P3  consensus length  ->  col2.txt   (get_fasta_length.sh inlined)
# ---------------------------------------------------------------------------

awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(substr(sequence_id, 2)" "sequence_length)}' cdhit.fa | grep "\S" > cdhit.fa.len
sort cdhit.fa.len | awk '{print $NF}' > col2.txt

echo ">>> [P3] DONE - consensus lengths calculated"


# ---------------------------------------------------------------------------
# P4  copy number in the genome  ->  col3.txt
# ---------------------------------------------------------------------------

if [ ! -f $genome.nin ]; then
    echo ">>> [P4] making blast database for the genome"
    makeblastdb -in $genome -dbtype nucl
fi

if [ ! -f genome.blast.o ]; then
    echo ">>> [P4] blasting library against the genome"
    blastn -query cdhit.fa -db $genome -num_threads $threads \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
        | awk '{OFS="\t"; if ($3 >= 80 && (($4/$13) > 0.5 )) {print $0,$4/$13}}' > genome.blast.o
fi

cat col1.txt genome.blast.o | awk '{print $1}' | sort | uniq -c | sort -k 2 | awk '{print $1-1}' > col3.txt

echo ">>> [P4] DONE - copy numbers estimated"


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

# Drop HSPs < MINLEN bp, rank remaining by BITSCORE, keep the best per family; record
# hit, pident, qcov, alignment length (bp), subject coverage (%), bitscore. High scov
# with low qcov -> a small RepBase element nested inside a large consensus.
awk -v m="$REPBASE_MINLEN" '$4>=m' repbase.blast \
  | sort -k1,1 -k7,7gr \
  | awk '!s[$1]++ {scov=($6>0?100*$4/$6:0); printf "%s\t%s\t%s\t%s\t%s\t%.1f\t%s\n",$1,$2,$3,$5,$4,scov,$7}' \
  | sort > repbase.best
join -a1 -e "0" -o 0,2.2,2.3,2.4,2.5,2.6,2.7 col1.txt repbase.best | tr ' ' '\t' | cut -f2-7 > col5.txt

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
# P6b  satellite screen  ->  col7.txt (sat_hit), col8.txt (sat_qcov)
#      Consensus = query, satellite lib = subject; blastn local + all HSPs, so a
#      consensus that is a tandem array of a monomer gets qcovs ~100%. -dust no
#      keeps low-complexity satellites from being masked away. High sat_qcov means
#      "this TE is (mostly) just a satellite".
# ---------------------------------------------------------------------------

if [ ! -f sat_db.nin ]; then
    echo ">>> [P6b] building blast db for the satellite library"
    makeblastdb -in $satlib -dbtype nucl -out sat_db
fi

if [ ! -f sat.blast ]; then
    echo ">>> [P6b] blasting library against satellites (dust off)"
    blastn -task blastn -dust no -query cdhit.fa -db sat_db -num_threads $threads \
        -outfmt "6 qseqid sseqid pident qcovs bitscore" -evalue 1e-5 > sat.blast
fi

# best satellite hit per family (highest bitscore) -> name + qcovs; join to col1 to stay aligned
sort -k1,1 -k5,5gr sat.blast | awk '!s[$1]++ {print $1"\t"$2"\t"$4}' | sort > sat.best
join -a1 -e "0" -o 0,2.2 col1.txt sat.best | tr ' ' '\t' | cut -f2 > col7.txt
join -a1 -e "0" -o 0,2.3 col1.txt sat.best | tr ' ' '\t' | cut -f2 > col8.txt

echo ">>> [P6b] DONE - `awk '$1!="0"' col7.txt | wc -l` families hit a satellite"


# ---------------------------------------------------------------------------
# P7  merge
# ---------------------------------------------------------------------------

echo -e "family\tlength\tcopies\tn_domains\trepbase_hit\trepbase_pident\trepbase_qcov\trepbase_alnlen\trepbase_scov\trepbase_bitscore\tdomains\tsat_hit\tsat_qcov" > final_priority.table.tab
paste -d "\t" col1.txt col2.txt col3.txt col4.txt col5.txt col6.txt col7.txt col8.txt >> final_priority.table.tab

echo ">>> [P7] DONE - final_priority.table.tab generated"
echo
echo "To list candidate families with no close RepBase match, ranked by copy number:"
echo "  awk -F'\t' 'NR>1 && (\$6==0 || \$6<70 || \$7<50)' final_priority.table.tab | sort -k3,3nr | less -S"
