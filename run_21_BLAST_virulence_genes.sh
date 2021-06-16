#!/bin/bash
###############################################################################
#21 Virulence gene blast

###############################################################################

echo "started.... Virulence gene blast ---------------------------------------"

###############################################################################
echo "provide list file (for e.g. all)"
ls list.*.txt | sed 's/ /\n/g'
read l
list=$(echo "list.$l.txt")
s=$(echo "_$l")
##-----------------------------------------------------------------------------
if [ ! -f /media/swapnil/share/databases/VFDB/VFDB_setB_pro.fas.phr ] ; then
echo "creating blast database"
makeblastdb -in /media/swapnil/share/databases/VFDB/VFDB_setB_pro.fas -parse_seqids -dbtype prot
else
echo "-------------------------------------------------------------------------------"
echo "VFDB Database versoin"
stat /media/swapnil/share/databases/VFDB/VFDB_setB_pro.fas.phr | awk 'NR >= 5 && NR <= 7'
echo "-------------------------------------------------------------------------------"
fi

(mkdir results/21_VFDB_blast_results"$s") > /dev/null 2>&1
(mkdir results/21_VFDB_blast_results"$s"/tmp) > /dev/null 2>&1

###############################################################################
for F1 in $(cat $list); do
    #(find results/ -name $F1.virulence-gene-blast.tmp -exec cp '{}' results/21_VFDB_blast_results"$s"/tmp/ \;) > /dev/null 2>&1
    #if [ ! -f results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.tmp ]; then
    echo "runnning... Virulence gene blast for $F1"
    blastp -db /media/swapnil/share/databases/VFDB/VFDB_setB_pro.fas -query results/08_annotation/raw_files/$F1/$F1.faa -out results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.tmp -max_target_seqs 1 -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
    #else
    #echo "Viruelence gene blast for $F1 already finished"
    #fi
    #--------------------
    # filtering the output for total-query-covergae >99% and total-protein identity >70%
    awk -F'\t' '{print ($8/$7*100) }' results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.tmp > results/21_VFDB_blast_results"$s"/tmp/$F1.total-query-covered.tmp
    paste results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.tmp results/21_VFDB_blast_results"$s"/tmp/$F1.total-query-covered.tmp > results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.total-query-covered.tmp
    awk -F'\t' '$14 > 70' results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.total-query-covered.tmp | awk -F'\t' '$16 > 95' > results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.total-query-covered.filtered.tmp
    #--------------------

    awk -f /home/swapnil/pipeline/tools/vlookup-VFDB.awk /media/swapnil/share/databases/VFDB/VFDB-annotations.txt results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.total-query-covered.filtered.tmp > results/21_VFDB_blast_results"$s"/tmp/$F1.tmp1.tmp

    sed -i 's/_/ /g' results/21_VFDB_blast_results"$s"/tmp/$F1.tmp1.tmp

    paste results/21_VFDB_blast_results"$s"/tmp/$F1.virulence-gene-blast.total-query-covered.filtered.tmp results/21_VFDB_blast_results"$s"/tmp/$F1.tmp1.tmp > results/21_VFDB_blast_results"$s"/$F1.VFDB-blast.csv

    ex -sc '1i|sseqid	qseqid	sstart	send	qstart	qend	slen	qlen	evalue	bitscore	length	mismatch	gaps	pident	qcovs	total-query-coverage	gene	protein	group	originated-from' -cx results/21_VFDB_blast_results"$s"/$F1.VFDB-blast.csv

    echo "finished... Virulence gene blast for $F1"

done

#------------------------------------------------------------------------------
#calcualte Abr-gene frequency

echo "summarising virulence gene blast analysis"

awk 'FNR>1' results/21_VFDB_blast_results"$s"/*.VFDB-blast.csv > results/21_VFDB_blast_results"$s"/all.VFDB-blast.tab
ex -sc '1i|sseqid	qseqid	sstart	send	qstart	qend	slen	qlen	evalue	bitscore	length	mismatch	gaps	pident	qcovs	total-query-coverage	gene	protein	group	originated-from' -cx results/21_VFDB_blast_results"$s"/all.VFDB-blast.tab
#unoconv -i FilterOptions=09,,system,1 -f xls results/21_VFDB_blast_results"$s"/all.VFDB-blast.tab

cp results/21_VFDB_blast_results"$s"/all.VFDB-blast.tab results/21_VFDB_blast_results"$s"/all.VFDB-blast.tab.tmp1
cp results/21_VFDB_blast_results"$s"/all.VFDB-blast.tab.tmp1 results/21_VFDB_blast_results"$s"/all.txt.1.tmp
sed -i 's/ /\t/g' results/21_VFDB_blast_results"$s"/all.VFDB-blast.tab

sed -i 's/ /_/g' results/21_VFDB_blast_results"$s"/all.txt.1.tmp
awk -F'\t' 'FNR >1 {print $17}' results/21_VFDB_blast_results"$s"/all.txt.1.tmp > results/21_VFDB_blast_results"$s"/all.txt.2.tmp
sed -i 's/  /===/' results/21_VFDB_blast_results"$s"/all.txt.2.tmp
sed -i 's/  /===/' results/21_VFDB_blast_results"$s"/all.txt.2.tmp
cat results/21_VFDB_blast_results"$s"/all.txt.2.tmp | sort | uniq > results/21_VFDB_blast_results"$s"/all.txt.3.tmp 

(rm results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.1.tmp) > /dev/null 2>&1 
for F2 in $(cat results/21_VFDB_blast_results"$s"/all.txt.3.tmp); do
V1=$(grep -c "$F2" results/21_VFDB_blast_results"$s"/all.txt.2.tmp)
echo $F2 $V1 >> results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.1.tmp
done

cat results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.1.tmp | sort -k 2,2rn results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.1.tmp > results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.2.tmp

awk -f /home/swapnil/pipeline/tools/vlookup-VFDB.2.awk /media/swapnil/share/databases/VFDB/VFDB-annotations.2.txt results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.2.tmp > results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.3.tmp

paste results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.2.tmp results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.3.tmp > results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv

ex -sc '1i|Virulence-gene frequency protein group originated-from' -cx results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv
sed -i 's/ /\t/g' results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv

sed -i 's/===/ /g' results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv
sed -i 's/_/ /g' results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv
sed -i 's/_/ /g' results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv
#unoconv -i FilterOptions=09,,system,1 -f xls results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv

echo "summarising virulence gene blast analysis completed"

###############################################################################
## Creating matrix out of VFDB blast results

echo "running matrix for VFDB"

for F1 in $(cat $list);do

awk -F'\t' 'FNR >1 {print $17}' results/21_VFDB_blast_results"$s"/$F1.VFDB-blast.csv | sed 's/ /_/g' > results/21_VFDB_blast_results"$s"/tmp/$F1.VFDB-blast.csv.tmp
	(rm results/21_VFDB_blast_results"$s"/$F1.gene-count.tmp1) > /dev/null 2>&1 
	for V1 in $(cat results/21_VFDB_blast_results"$s"/all.txt.3.tmp);do
	V2=$(awk '{count[$1]++} END {print count["'$V1'"]}' results/21_VFDB_blast_results"$s"/tmp/$F1.VFDB-blast.csv.tmp)
	echo $V2 >> results/21_VFDB_blast_results"$s"/$F1.gene-count.tmp1
	done
awk '{for (i=1; i<= NF; i++) {if($i > 1) { $i=1; } } print }' results/21_VFDB_blast_results"$s"/$F1.gene-count.tmp1 > results/21_VFDB_blast_results"$s"/$F1.gene-count.tmp2
ex -sc '1i|'$F1'' -cx results/21_VFDB_blast_results"$s"/$F1.gene-count.tmp1
ex -sc '1i|'$F1'' -cx results/21_VFDB_blast_results"$s"/$F1.gene-count.tmp2
sed -i -e 's/^$/0/' results/21_VFDB_blast_results"$s"/$F1.gene-count.tmp1
sed -i -e 's/^$/0/' results/21_VFDB_blast_results"$s"/$F1.gene-count.tmp2
done

cp results/21_VFDB_blast_results"$s"/all.txt.3.tmp results/21_VFDB_blast_results"$s"/all.txt.4.tmp
sed -i 's/$/:c/' results/21_VFDB_blast_results"$s"/all.txt.4.tmp
sed -i '1i\===\' results/21_VFDB_blast_results"$s"/all.txt.4.tmp
paste results/21_VFDB_blast_results"$s"/all.txt.4.tmp results/21_VFDB_blast_results"$s"/*.gene-count.tmp1 > results/21_VFDB_blast_results"$s"/matrix.csv.tmp1
paste results/21_VFDB_blast_results"$s"/all.txt.4.tmp results/21_VFDB_blast_results"$s"/*.gene-count.tmp2 > results/21_VFDB_blast_results"$s"/matrix_1-0.csv.tmp2

##-----------------------------------------------------------------------------
## transpose
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' results/21_VFDB_blast_results"$s"/matrix.csv.tmp1 > results/21_VFDB_blast_results"$s"/matrix.csv

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' results/21_VFDB_blast_results"$s"/matrix_1-0.csv.tmp2 > results/21_VFDB_blast_results"$s"/matrix_1-0.csv

##-----------------------------------------------------------------------------

sed -i 's/===//g' results/21_VFDB_blast_results"$s"/matrix.csv
sed -i 's/===//g' results/21_VFDB_blast_results"$s"/matrix_1-0.csv
sed -i 's/ /,/g' results/21_VFDB_blast_results"$s"/matrix.csv
sed -i 's/ /,/g' results/21_VFDB_blast_results"$s"/matrix_1-0.csv

rm results/21_VFDB_blast_results"$s"/*.tmp1
rm results/21_VFDB_blast_results"$s"/*.tmp2
rm results/21_VFDB_blast_results"$s"/Virulence_gene_frequency.csv.1.tmp

echo "finished matrix for VFDB"

###############################################################################

echo "Completed.. Virulence gene blast ---------------------------------------"

###############################################################################
