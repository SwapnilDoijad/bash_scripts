!/bin/bash
###############################################################################
#6 16SrRNA_gene

###############################################################################

echo "started...... 16SrRNA gene ---------------------------------------------"

###############################################################################
    #makeblastdb -in /media/swapnil/share/databases/16SMicrobial/16SMicrobial.fasta -parse_seqids -dbtype nucl

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/06_16SrDNA) > /dev/null 2>&1
    (mkdir results/06_16SrDNA/tmp) > /dev/null 2>&1 
    (mkdir results/06_16SrDNA/raw_files/) > /dev/null 2>&1 
    (mkdir results/06_16SrDNA/raw_files/tmp) > /dev/null 2>&1 
###############################################################################

for F1 in $(cat $list); do

echo "running 16S rDNA blast for $F1"

#perl /home/swapnil/tools/rnammer/rnammer -S bac -m ssu -f results/06_16SrDNA/raw_files/$F1.16SrDNA.fasta results/04_assembly/all_fasta/$F1.fasta 

(/home/swapnil/tools/barrnap/bin/barrnap --quiet --outseq results/06_16SrDNA/raw_files/$F1.rDNA.out.fasta results/04_assembly/all_fasta/$F1.fasta ) > /dev/null 2>&1 

/home/swapnil/pipeline/tools/fastagrep.pl -X 16S_rRNA results/06_16SrDNA/raw_files/$F1.rDNA.out.fasta > results/06_16SrDNA/raw_files/$F1.16SrDNA.fasta

blastn -db /media/swapnil/share/databases/16SMicrobial/16SMicrobial.fasta -query results/06_16SrDNA/raw_files/$F1.16SrDNA.fasta -out results/06_16SrDNA/raw_files/tmp/results.tab -max_target_seqs 5 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs qseq"

awk -f /home/swapnil/pipeline/tools/vlookup-16S.awk /media/swapnil/share/databases/16SMicrobial/header.tab results/06_16SrDNA/raw_files/tmp/results.tab > results/06_16SrDNA/raw_files/tmp/headers.tmp

paste results/06_16SrDNA/raw_files/tmp/results.tab results/06_16SrDNA/raw_files/tmp/headers.tmp > results/06_16SrDNA/raw_files/tmp/$F1.16S_rDNA_BLAST_results.csv.tmp

echo "finished 16S rDNA blast for $F1"

done

cat results/06_16SrDNA/raw_files/tmp/*.16S_rDNA_BLAST_results.csv.tmp > results/06_16SrDNA/tmp/16S_rDNA_BLAST_results.csv

ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs qseq 16S' -cx results/06_16SrDNA/tmp/16S_rDNA_BLAST_results.csv

sed -i "1 s/ /\t/g" results/06_16SrDNA/tmp/16S_rDNA_BLAST_results.csv

unoconv -i FilterOptions=09,,system,1 -f xls -o results/06_16SrDNA/16S_rDNA_BLAST_results results/06_16SrDNA/tmp/16S_rDNA_BLAST_results.csv

###############################################################################

echo "finished..... 16SrRNA gene ---------------------------------------------"

###############################################################################
