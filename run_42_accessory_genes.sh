#!/bin/bash
###############################################################################
if [ ! -d results/42_accessory_gene_analysis ]; then
mkdir results/42_accessory_gene_analysis
fi

if [ ! -d results/42_accessory_gene_analysis/tmp ]; then
mkdir results/42_accessory_gene_analysis/tmp
fi

if [ ! -d results/42_accessory_gene_analysis/tmp/prophage ]; then
mkdir results/42_accessory_gene_analysis/tmp/prophage
fi

cp results/17_roary/roary_results/gene_presence_absence.csv results/42_accessory_gene_analysis/tmp/

#------------------------------------------------------------------------------

## count number of isolates

sed '/^\s*$/d' list.txt | wc -l  | awk '{print $1}' > results/42_accessory_gene_analysis/tmp/number_of_isolates.txt

V1=$(cat results/42_accessory_gene_analysis/tmp/number_of_isolates.txt) ;
V2=$(( $V1 - 1 )) ;

awk -vFPAT='([^,]*)|("[^"]+")' -vOFS=, ' $4 <= '$V2' {print $0} ' results/42_accessory_gene_analysis/tmp/gene_presence_absence.csv > results/42_accessory_gene_analysis/tmp/accessory.tmp

awk -vFPAT='([^,]*)|("[^"]+")' -vOFS=, ' FNR ==1 {print $0}' results/17_roary/roary_results/gene_presence_absence.csv > results/42_accessory_gene_analysis/tmp/header.txt

cat results/42_accessory_gene_analysis/tmp/header.txt results/42_accessory_gene_analysis/tmp/accessory.tmp > results/42_accessory_gene_analysis/tmp/accessory.csv

awk -vFPAT='([^,]*)|("[^"]+")' -vOFS=, '{ s = ""; for (i = 15; i <= NF; i++) s = s $i " "; print s }' results/42_accessory_gene_analysis/tmp/accessory.tmp > results/42_accessory_gene_analysis/tmp/accessory.only_genes.tmp

awk '{print $1}' results/42_accessory_gene_analysis/tmp/accessory.only_genes.tmp > results/42_accessory_gene_analysis/tmp/accessory.only_genes.one_column.tmp


#------------------------------------------------------------------------------
## prophage

echo "identifying prophage genes"

for F1 in $(cat list.txt); do

echo "Identifying prophage genes in $F1"

cat results/30_prophage_finder/results/$F1/$F1/$F1.blastt-best-matches-ids >> results/42_accessory_gene_analysis/tmp/prophage/all.blastt-best-matches-ids ;
done

sed -i 's/\.t01//g' results/42_accessory_gene_analysis/tmp/prophage/all.blastt-best-matches-ids ;


for F2 in $(cat results/42_accessory_gene_analysis/tmp/prophage/all.blastt-best-matches-ids) ; do
echo "writing prophage genes for $F1"
sed -i 's/'$F2'/'$F2'_prophage/g' results/42_accessory_gene_analysis/tmp/accessory.csv
done


#------------------------------------------------------------------------------


## BLAST gene to detect Plasmid gene 

echo "identifying plasmid genes"

for F1 in $(cat list.txt); do

cat results/08_annotation/$F1/*.ffn >> results/42_accessory_gene_analysis/tmp/all.ffn

done
echo "extracting plasmid genes"

perl  /home/swapnil/pipeline/tools/fastagrep.pl -f results/42_accessory_gene_analysis/tmp/accessory.only_genes.one_column.tmp results/42_accessory_gene_analysis/tmp/all.ffn > results/42_accessory_gene_analysis/tmp/accessory.only_genes.one_column.tmp.ffn

sed -i 's/ .*//g' results/42_accessory_gene_analysis/tmp/accessory.only_genes.one_column.tmp.ffn

echo "BLASTing plasmid genes"

blastn -db /home/swapnil/pipeline/tools/databases/plasmid/refseq_plasmid.fasta -query results/42_accessory_gene_analysis/tmp/accessory.only_genes.one_column.tmp.ffn -max_target_seqs 1 -max_hsps 1 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" > results/42_accessory_gene_analysis/tmp/blast.plasmid.tmp

awk '$15 >= 70 {print $2}' results/42_accessory_gene_analysis/tmp/blast.plasmid.tmp > results/42_accessory_gene_analysis/tmp/blast.plasmid.genes.70filter.tmp


for F3 in $(cat results/42_accessory_gene_analysis/tmp/blast.plasmid.genes.70filter.tmp) ; do
echo "writing plasmid genes for $F3"
sed -i 's/'$F3'/'$F3'_plasmid/g' results/42_accessory_gene_analysis/tmp/accessory.csv
done


#------------------------------------------------------------------------------

## count total, highest and lowest number of accessory genes

echo "counting accessory gene numbers"

V3=$(awk '{print NR}' results/42_accessory_gene_analysis/tmp/accessory.csv | tail -1)

V4=$(( $V3 - 1 ))

echo total_number_of_accessory_genes $V3 > results/42_accessory_gene_analysis/number_of_accessory_gens.csv
echo "  " >> results/42_accessory_gene_analysis/number_of_accessory_gens.csv
echo "Isolate accessory_genes" >> results/42_accessory_gene_analysis/number_of_accessory_gens.csv

for F1 in $(cat list.txt); do

awk -F',' 'NR==1{for (i=1;i<=NF;i++) a[i]=$i; next} {for (i=1;i<=NF;i++) {print $i > a[i]".tmptmp"}}' results/42_accessory_gene_analysis/tmp/accessory.csv ;

sed -i '/^\s*$/d' $F1.tmptmp

V5=$(wc -l $F1.tmptmp | awk '{print $1}') ;
V6=$(( $V5 - 1 )) ;

echo $F1 $V6 >> results/42_accessory_gene_analysis/number_of_accessory_gens.csv

rm *.tmptmp

done

unoconv -e FilterOptions=44,34,76 -f xls -o results/42_accessory_gene_analysis/accessory.xls results/42_accessory_gene_analysis/tmp/accessory.csv




