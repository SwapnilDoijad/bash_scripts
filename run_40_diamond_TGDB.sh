#!/bin/bash
###############################################################################
#40 transporter gene blast

###############################################################################

echo "started.... transporter gene blast --------------------------------------"

###############################################################################

echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

echo "need to add suffix for the result folder? type 'y' or press enter"
read suffix_answer
if [ "$suffix_answer" == "y" ]; then
s=$(echo "_$l")
fi

##-----------------------------------------------------------------------------

if [ ! -f /home/swapnil/pipeline/tools/databases/tcdb/nr.dmnd ] ; then
diamond makedb --in /home/swapnil/pipeline/tools/databases/tcdb/030519_tcdb.faa -d /home/swapnil/pipeline/tools/databases/tcdb/nr
fi

if [ ! -d results/40_tcdb_diamond_results"$s" ]; then
mkdir results/40_tcdb_diamond_results"$s"
fi

if [ ! -d results/40_tcdb_diamond_results"$s"/tmp ]; then
mkdir results/40_tcdb_diamond_results"$s"/tmp
fi

(rm results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.1.tmp)> /dev/null 2>&1

#--------------------
for F1 in $(cat $list); do
echo "runnning... transporter gene blast for $F1"
#--------------------

## check if blast is already carried out
#(find results/ -name $F1.tcdb-blast.tmp -exec cp '{}' results/21_tcdb_blast_results"$s"/ \;) > /dev/null 2>&1

#if [ ! -f results/21_tcdb_blast_results/$F1.tcdb-blast.tmp ]; then

(diamond blastp -b 8 -d /home/swapnil/pipeline/tools/databases/tcdb/nr -q results/08_annotation/raw_files/$F1/$F1.faa --id 50 --query-cover 50 --max-target-seqs 1 -o results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.tmp -f 6 sseqid qseqid slen qlen evalue bitscore length pident nident mismatch gaps)> /dev/null 2>&1 

#--------------------
# filtering the output for total-query-covergae >99% and total-protein identity >70%
awk -F'\t' '{print ($7/$4*100) }' results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.tmp > results/40_tcdb_diamond_results"$s"/tmp/$F1.total-query-covered.tmp
paste results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.tmp results/40_tcdb_diamond_results"$s"/tmp/$F1.total-query-covered.tmp > results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.total-query-covered.filtered.tmp
#--------------------

awk -F'|' '{print $3}' results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.total-query-covered.filtered.tmp > results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.total-query-covered.filtered.2.tmp

awk -f /home/swapnil/pipeline/tools/vlookup-tcdb.awk /home/swapnil/pipeline/tools/databases/tcdb/tcdb-annotations.txt results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.total-query-covered.filtered.2.tmp > results/40_tcdb_diamond_results"$s"/tmp/$F1.tmp1.tmp

sed -i 's/_/ /g' results/40_tcdb_diamond_results"$s"/tmp/$F1.tmp1.tmp

paste results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.total-query-covered.filtered.tmp results/40_tcdb_diamond_results"$s"/tmp/$F1.tmp1.tmp > results/40_tcdb_diamond_results"$s"/$F1.tcdb-blast.csv

ex -sc '1i|sseqid	qseqid	slen	qlen	evalue	bitscore	length	pident	nident	mismatch	gaps	total_QrCov	name' -cx results/40_tcdb_diamond_results"$s"/$F1.tcdb-blast.csv

ssconvert results/40_tcdb_diamond_results"$s"/$F1.tcdb-blast.csv results/40_tcdb_diamond_results"$s"/$F1.tcdb-blast.xlsx
#fi
echo "finished... tcdb gene blast for $F1"
done

#------------------------------------------------------------------------------
#calcualte gene frequency

#echo "summarising tcdb gene blast analysis"

#awk 'FNR>1' results/40_tcdb_diamond_results"$s"/*.tcdb-blast.csv > results/40_tcdb_diamond_results"$s"/all.tcdb-blast.tab
#ex -sc '1i|sseqid	qseqid	sstart	send	qstart	qend	slen	qlen	evalue	bitscore	length	mismatch	gaps	pident	qcovs	total-query-coverage	gene	protein	group	originated-from' -cx results/40_tcdb_diamond_results"$s"/all.tcdb-blast.tab
#unoconv -i FilterOptions=09,,system,1 -f xls results/40_tcdb_diamond_results"$s"/all.tcdb-blast.tab

#cp results/40_tcdb_diamond_results"$s"/all.tcdb-blast.tab results/40_tcdb_diamond_results"$s"/all.tcdb-blast.tab.tmp1
#cp results/40_tcdb_diamond_results"$s"/all.tcdb-blast.tab.tmp1 results/40_tcdb_diamond_results"$s"/all.txt.1.tmp
#sed -i 's/ /\t/g' results/40_tcdb_diamond_results"$s"/all.tcdb-blast.tab

#sed -i 's/ /_/g' results/40_tcdb_diamond_results"$s"/all.txt.1.tmp
#awk -F'\t' 'FNR >1 {print $13}' results/40_tcdb_diamond_results"$s"/all.txt.1.tmp > results/40_tcdb_diamond_results"$s"/all.txt.2.tmp
#sed -i 's/  /===/' results/40_tcdb_diamond_results"$s"/all.txt.2.tmp
#sed -i 's/  /===/' results/40_tcdb_diamond_results"$s"/all.txt.2.tmp
#cat results/40_tcdb_diamond_results"$s"/all.txt.2.tmp | sort | uniq > results/40_tcdb_diamond_results"$s"/all.txt.3.tmp 

#for F2 in $(cat results/40_tcdb_diamond_results"$s"/all.txt.3.tmp); do
#V1=$(grep -c "$F2" results/40_tcdb_diamond_results"$s"/all.txt.2.tmp)
#echo $F2 $V1 >> results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.1.tmp
#done

#cat results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.1.tmp | sort -k 2,2rn results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.1.tmp > results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.2.tmp

#awk -f /home/swapnil/pipeline/tools/vlookup-tcdb.2.awk /home/swapnil/tools/databases/tcdb/tcdb-annotations.2.txt results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.2.tmp > results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.3.tmp

#paste results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.2.tmp results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.3.tmp > results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv

#ex -sc '1i|tcdb frequency protein group originated-from' -cx results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv
#sed -i 's/ /\t/g' results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv

#sed -i 's/===/ /g' results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv
#sed -i 's/_/ /g' results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv
#sed -i 's/_/ /g' results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv
#unoconv -i FilterOptions=09,,system,1 -f xls results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv

#echo "summarising tcdb gene blast analysis completed"

###############################################################################
## Creating matrix out of tcdb blast results

#echo "running matrix for tcdb"

#for F1 in $(cat $list);do

#awk -F'\t' 'FNR >1 {print $13}' results/40_tcdb_diamond_results"$s"/$F1.tcdb-blast.csv | sed 's/ /_/g' > results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.csv.tmp
#	for V1 in $(cat results/40_tcdb_diamond_results"$s"/all.txt.3.tmp);do
#	V2=$(awk '{count[$1]++} END {print count["'$V1'"]}' results/40_tcdb_diamond_results"$s"/tmp/$F1.tcdb-blast.csv.tmp)
#	echo $V2 >> results/40_tcdb_diamond_results"$s"/$F1.gene-count.tmp1
#	done
#awk '{for (i=1; i<= NF; i++) {if($i > 1) { $i=1; } } print }' results/40_tcdb_diamond_results"$s"/$F1.gene-count.tmp1 > results/40_tcdb_diamond_results"$s"/$F1.gene-count.tmp2
#ex -sc '1i|'$F1'' -cx results/40_tcdb_diamond_results"$s"/$F1.gene-count.tmp1
#ex -sc '1i|'$F1'' -cx results/40_tcdb_diamond_results"$s"/$F1.gene-count.tmp2
#sed -i -e 's/^$/0/' results/40_tcdb_diamond_results"$s"/$F1.gene-count.tmp1
#sed -i -e 's/^$/0/' results/40_tcdb_diamond_results"$s"/$F1.gene-count.tmp2
#done

#cp results/40_tcdb_diamond_results"$s"/all.txt.3.tmp results/40_tcdb_diamond_results"$s"/all.txt.4.tmp
#sed -i 's/$/:c/' results/40_tcdb_diamond_results"$s"/all.txt.4.tmp
#sed -i '1i\===\' results/40_tcdb_diamond_results"$s"/all.txt.4.tmp
#paste results/40_tcdb_diamond_results"$s"/all.txt.4.tmp results/40_tcdb_diamond_results"$s"/*.gene-count.tmp1 > results/40_tcdb_diamond_results"$s"/matrix.csv.tmp1
#paste results/40_tcdb_diamond_results"$s"/all.txt.4.tmp results/40_tcdb_diamond_results"$s"/*.gene-count.tmp2 > results/40_tcdb_diamond_results"$s"/matrix_1-0.csv.tmp2

##-----------------------------------------------------------------------------
## transpose
#awk '
#{ 
#    for (i=1; i<=NF; i++)  {
#        a[NR,i] = $i
#    }
#}
#NF>p { p = NF }
#END {    
#    for(j=1; j<=p; j++) {
#        str=a[1,j]
#        for(i=2; i<=NR; i++){
#            str=str" "a[i,j];
#        }
#        print str
#    }
#}' results/40_tcdb_diamond_results"$s"/matrix.csv.tmp1 > results/40_tcdb_diamond_results"$s"/matrix.csv

#awk '
#{ 
#    for (i=1; i<=NF; i++)  {
#        a[NR,i] = $i
#    }
#}
#NF>p { p = NF }
#END {    
#    for(j=1; j<=p; j++) {
#        str=a[1,j]
#        for(i=2; i<=NR; i++){
#            str=str" "a[i,j];
#        }
#        print str
#    }
#}' results/40_tcdb_diamond_results"$s"/matrix_1-0.csv.tmp2 > results/40_tcdb_diamond_results"$s"/matrix_1-0.csv
##-----------------------------------------------------------------------------

#sed -i 's/===//g' results/40_tcdb_diamond_results"$s"/matrix.csv
#sed -i 's/===//g' results/40_tcdb_diamond_results"$s"/matrix_1-0.csv
#sed -i 's/ /,/g' results/40_tcdb_diamond_results"$s"/matrix.csv
#sed -i 's/ /,/g' results/40_tcdb_diamond_results"$s"/matrix_1-0.csv

#rm results/40_tcdb_diamond_results"$s"/*.tmp1
#rm results/40_tcdb_diamond_results"$s"/*.tmp2
#rm results/40_tcdb_diamond_results"$s"/tcdb_gene_frequency.csv.1.tmp

#echo "finished matrix for tcdb"

###############################################################################

echo "Completed.. tcdb blast --------------------------------------------------"

###############################################################################
