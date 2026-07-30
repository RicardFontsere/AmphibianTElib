#!/bin/bash

set -e

if [ $# -ne 5 ]
then
    echo -e "\nusage: `basename $0` <TEtrimmer_lib.fa> <genome.fa> <pfam_db_dir> <repbase_blastdb> <threads>\n"
    echo -e "DESCRIPTION:\tRuns a pipeline that: 1) reduces sequence redundancy from the TEtrimmer library"
    echo -e "\t\tusing cd-hit-est; 2) extracts family names; 3) calculates consensus length; 4) makes a rough"
    echo -e "\t\testimate of genome copy number; 5) screens each family against RepBase for novelty;"
    echo -e "\t\t6) predicts ORFs and identifies Pfam domains; 7) merges into final_priority.table.tab.\n"
    echo -e "REQUIRES:\tcd-hit, blast, EMBOSS (getorf), pfam_scan.pl and the Pfam database.\n"
    echo -e "INPUT:\t\t<TEtrimmer_lib.fa>\tTEtrimmer consensus library (e.g. TEtrimmer_consensus_merged.fasta)"
    echo -e "\t\t<genome.fa>\t\tgenome used to predict the library"
    echo -e "\t\t<pfam_db_DIR>\t\tpath to the Pfam database directory"
    echo -e "\t\t<repbase_blastdb>\tpath to a PRE-BUILT RepBase blast nucl database (db prefix, already makeblastdb'd)"
    echo -e "\t\t<threads>\t\tnumber of CPU threads"
    echo -e "OUTPUT:\t\tfinal_priority.table.tab, with columns:"
    echo -e "\t\tfamily | length | copies | n_domains | repbase_hit | repbase_pident | repbase_qcov | domains\n"
    exit
fi

rmout=$1
genome=$2
pfamdb=$3
repbase=$4
threads=$5


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

if [ ! -f repbase.blast ]; then
    echo ">>> [P5] blasting library against pre-built RepBase db"
    blastn -task dc-megablast -query cdhit.fa -db $repbase -num_threads $threads \
        -outfmt "6 qseqid sseqid pident qcovs bitscore" -evalue 1e-5 > repbase.blast
fi

sort -k1,1 -k5,5gr repbase.blast | awk '!s[$1]++ {print $1"\t"$2"\t"$3"\t"$4}' | sort > repbase.best
join -a1 -e "-" -o 0,2.2,2.3,2.4 col1.txt repbase.best | tr ' ' '\t' | cut -f2,3,4 > col5.txt

echo ">>> [P5] DONE - `awk '$1!="-"' col5.txt | wc -l` families have a RepBase hit"


# ---------------------------------------------------------------------------
# P6  ORF prediction and Pfam domains  ->  col4.txt (counts), col6.txt (names)
# ---------------------------------------------------------------------------

if [ ! -f cdhit.orf ]; then
    echo ">>> [P6] predicting ORFs"
    getorf -sequence cdhit.fa -outseq cdhit.orf -minsize 300
fi

if [ ! -f pfam.results ]; 
    then
    export PERL5LIB=/user/brussel/109/vsc10945/home/scratch/Software/PfamScan:$PERL5LIB
    echo ">>> [P6] running pfam_scan.pl, this can take some time"
    /user/brussel/109/vsc10945/home/scratch/Software/PfamScan/pfam_scan.pl -fasta cdhit.orf -dir $pfamdb -cpu $threads > pfam.results
fi

# domain counts
awk '{if ($6~/^PF/) {print $1}}' < pfam.results | sed 's/_/\//2;s/_/ /2' | awk '{print $1}' | sort > pf.domains.count
cat col1.txt pf.domains.count | sort | uniq -c | sort -k 2 | awk '{print $1-1}' > col4.txt

# domain names, one comma-separated list per family
awk '$6~/^PF/ {n=$1; sub(/_[0-9]+$/,"",n); k=n"\t"$7; if(!s[k]++) d[n]=(n in d ? d[n]","$7 : $7)}
     END {for(f in d) print f"\t"d[f]}' pfam.results | sort > pf.domains.names
join -a1 -e "none" -o 0,2.2 col1.txt pf.domains.names | tr ' ' '\t' | cut -f2 > col6.txt

echo ">>> [P6] DONE - Pfam searches finished"


# ---------------------------------------------------------------------------
# P7  merge
# ---------------------------------------------------------------------------

echo -e "family\tlength\tcopies\tn_domains\trepbase_hit\trepbase_pident\trepbase_qcov\tdomains" > final_priority.table.tab
paste -d "\t" col1.txt col2.txt col3.txt col4.txt col5.txt col6.txt >> final_priority.table.tab

echo ">>> [P7] DONE - final_priority.table.tab generated"
echo
echo "To list candidate families with no close RepBase match, ranked by copy number:"
echo "  awk -F'\t' 'NR>1 && (\$6==\"-\" || \$6<70 || \$7<50)' final_priority.table.tab | sort -k3,3nr | less -S"
