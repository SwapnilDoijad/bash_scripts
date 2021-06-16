#!/bin/bash
###############################################################################
#23 BLAST custom gene set
###############################################################################
echo "started.... BLAST custom gene set --------------------------------------"
###############################################################################

echo "provide list file (for e.g. all)"
echo "-------------------------------------------------------------------------"
ls list.*.txt | sed 's/ /\n/g'
echo "-------------------------------------------------------------------------"
read l
list=$(echo "list.$l.txt")

echo "provide gene files (for e.g. data/genes/bcnA.fasta)"
echo "-------------------------------------------------------------------------"
ls data/genes/*.fasta | sed 's/ /\n/g'
echo "-------------------------------------------------------------------------"
read ref_genes
gene=$(echo $ref_genes | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}')

(mkdir results/23_Alignment_$gene) > /dev/null 2>&1
(mkdir results/23_Alignment_$gene/tmp) > /dev/null 2>&1
(mkdir data/genes) > /dev/null 2>&1

echo "-------------------------------------------------------------------------"
makeblastdb -in $ref_genes -parse_seqids -dbtype nucl
echo "-------------------------------------------------------------------------"

echo "for DNA type d or Protein type p"
read answer
###############################################################################

if [ "$answer" = "d" ]; then

    for F1 in $(cat $list); do
    if [ ! -f results/23_Alignment_$gene/tmp/$F1.custom-gene-blast.tmp ]; then
        echo "runnning... BLAST custom gene set $F1"
        blastn -db $ref_genes -query results/04_assembly/all_fasta/$F1.fasta -out results/23_Alignment_$gene/tmp/$F1.custom-gene-blast.tmp -max_target_seqs 100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
        #blastn -db $ref_genes -query results/04_assembly/all_fasta/$F1.fasta -out results/23_Alignment_$gene/tmp/$F1.custom-gene-blast.tmp -max_target_seqs 1 -evalue 1e-100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
        ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs' -cx results/23_Alignment_$gene/tmp/$F1.custom-gene-blast.tmp
        sed -i '1s/ /\t/g' results/23_Alignment_$gene/tmp/$F1.custom-gene-blast.tmp
        else
        echo "$F1 is already finished"
        fi
    done

elif [ "$answer" = "p" ]; then
    
    for F1 in $(cat $list); do
    if [ ! -f results/23_Alignment_$gene/tmp/$F1.$gene.total-query-covered.filtered.tmp ]; then
    echo "runnning... DIAMOND for custom gene set $F1"
    diamond makedb --in $ref_genes -d data/genes
    (diamond blastp -b 8 -d data/genes -q results/08_annotation/raw_files/$F1/$F1.faa --id 80 --query-cover 70 --max-target-seqs 1 -o results/23_Alignment_$gene/tmp/$F1.$gene.tmp -f 6 sseqid qseqid slen qlen evalue bitscore length pident nident mismatch gaps )> /dev/null 2>&1 
    #--------------------
    # filtering the output for total-query-covergae >99% and total-protein identity >70%
    awk -F'\t' '{print ($7/$4*100) }' results/23_Alignment_$gene/tmp/$F1.$gene.tmp > results/23_Alignment_$gene/tmp/$F1.$gene.total-query-covered.tmp
    paste results/23_Alignment_$gene/tmp/$F1.$gene.tmp results/23_Alignment_$gene/tmp/$F1.$gene.total-query-covered.tmp > results/23_Alignment_$gene/tmp/$F1.$gene.total-query-covered.filtered.tmp
    ex -sc '1i|sseqid	qseqid	slen	qlen	evalue	bitscore	length	pident	nident	mismatch	gaps    total_QC' -cx results/23_Alignment_$gene/tmp/$F1.$gene.total-query-covered.filtered.tmp
    else
    echo "$F1 is already finished"
    fi
    done

else
    echo "Unknown parameter"
fi

    awk 'FNR>1' results/23_Alignment_$gene/tmp/*.tmp > results/23_Alignment_$gene/all.csv
    ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs' -cx results/23_Alignment_$gene/all.csv
    sed -i '1s/ /\t/g' results/23_Alignment_$gene/all.csv
###############################################################################
echo "Completed.... BLAST custom gene set ------------------------------------"
###############################################################################
