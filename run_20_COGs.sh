#!/bin/bash
###############################################################################
## step-20 COGs

###############################################################################
if [ ! -d results/20_COGs ]; then
mkdir results/20_COGs
fi

for F1 in $(cat list.txt); do

echo "running COG analysis for $F1"

if [ ! -d results/20_COGs/$F1 ]; then
mkdir results/20_COGs/$F1
fi

mkdir results/20_COGs/$F1/BLASTff
mkdir results/20_COGs/$F1/BLASTno
mkdir results/20_COGs/$F1/BLASTss
mkdir results/20_COGs/$F1/BLASTcogn
mkdir results/20_COGs/$F1/tmp

(makeblastdb -in results/08_annotation/all_faa/$F1.faa -dbtype prot -out results/20_COGs/$F1/tmp/$F1) > /dev/null 2>&1

(psiblast -query results/08_annotation/all_faa/$F1.faa -db results/20_COGs/$F1/tmp/$F1 -show_gis -outfmt 7 -num_descriptions 10 -num_alignments 10 -dbsize 100000000 -comp_based_stats F -seg no -out results/20_COGs/$F1/BLASTss/QuerySelf.tab) > /dev/null 2>&1

(psiblast -query results/08_annotation/all_faa/$F1.faa -db /home/swapnil/pipeline/tools/COGs/COGs -show_gis -outfmt 7 -num_descriptions 1000 -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out results/20_COGs/$F1/BLASTno/QueryCOGs.tab) > /dev/null 2>&1

(psiblast -query results/08_annotation/all_faa/$F1.faa -db /home/swapnil/pipeline/tools/COGs/COGs -show_gis -outfmt 7 -num_descriptions 1000 -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out results/20_COGs/$F1/BLASTff/QueryCOGs.tab) > /dev/null 2>&1

grep ">" results/08_annotation/all_faa/$F1.faa | awk '{print $1",""'$F1'"}' | sed 's/>//g' > results/20_COGs/$F1/tmp/$F1.gene_list.csv

cat results/20_COGs/$F1/tmp/$F1.gene_list.csv /home/swapnil/pipeline/tools/COGs/COGs.csv > results/20_COGs/$F1/tmp/$F1.gene_list.COGs_list.tmp.csv

COGmakehash -i=results/20_COGs/$F1/tmp/$F1.gene_list.COGs_list.tmp.csv -o=results/20_COGs/$F1/BLASTcogn -s="," -n=1 

COGreadblast -d=results/20_COGs/$F1/BLASTcogn -u=results/20_COGs/$F1/BLASTno -f=results/20_COGs/$F1/BLASTff -s=results/20_COGs/$F1/BLASTss -e=0.1 -q=2 -t=2

COGcognitor -i=results/20_COGs/$F1/BLASTcogn -t=/home/swapnil/pipeline/tools/COGs/COGs.csv -q=results/20_COGs/$F1/tmp/$F1.gene_list.csv -o=results/20_COGs/$F1/tmp/$F1.COGs.csv.tmp

grep ",COG" results/20_COGs/$F1/tmp/$F1.COGs.csv.tmp | awk -F',' '{print $1, $6}' > results/20_COGs/$F1/tmp/$F1.COGs.csv

awk '{print $2}' results/20_COGs/$F1/tmp/$F1.COGs.csv > results/20_COGs/$F1/tmp/$F1.COGs.COGs_list.tmp

awk -f /home/swapnil/pipeline/tools/vlookup-COGs.awk /home/swapnil/pipeline/tools/COGs/cognames2003-2014.tab results/20_COGs/$F1/tmp/$F1.COGs.COGs_list.tmp > results/20_COGs/$F1/tmp/$F1.COGs.COGs_list.2.tmp
sed -i '/^\s*$/d' results/20_COGs/$F1/tmp/$F1.COGs.COGs_list.2.tmp

	for F2 in $(cat /home/swapnil/pipeline/tools/COGs/functional_category_list.tab);do
	V1=$(grep -c "$F2" results/20_COGs/$F1/tmp/$F1.COGs.COGs_list.2.tmp)
	echo $V1 >> results/20_COGs/$F1/tmp/$F1.functional_category_list.count.tmp
	done
paste /home/swapnil/pipeline/tools/COGs/functional_category_list.tab /home/swapnil/pipeline/tools/COGs/functional_category_names.tab results/20_COGs/$F1/tmp/$F1.functional_category_list.count.tmp > results/20_COGs/$F1/tmp/$F1.functional_category.csv
sed -i 's/ /\t/g' results/20_COGs/$F1/tmp/$F1.functional_category.csv
sed -i 's/_/ /g' results/20_COGs/$F1/tmp/$F1.functional_category.csv
unoconv -i FilterOptions=09,,system,1 -f xls results/20_COGs/$F1/tmp/$F1.functional_category.csv

awk -F'\t' -v OFS='\t' '{print $2, $3}' results/20_COGs/$F1/tmp/$F1.functional_category.csv > results/20_COGs/$F1/tmp/$F1.functional_category.for_gnuplot.csv.tmp

awk -F'\t' '{sub($1, "\"&\""); print}' results/20_COGs/$F1/tmp/$F1.functional_category.for_gnuplot.csv.tmp > results/20_COGs/functional_category.for_gnuplot.dat

gnuplot /home/swapnil/pipeline/tools/COGs/simple.barchart.gnuplot 

mv barchart.png results/20_COGs/$F1/

rm results/20_COGs/functional_category.for_gnuplot.dat
rm cognitor.log
rm conflict.txt

echo "finished COG analysis for $F1"

done
