#22 CARD-ABR blast

###############################################################################

echo "started.... CARD-ABR blast ---------------------------------------------"

###############################################################################

echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

echo "want to suffix files? type "y" else press enter to continue"
read answer
if [ "$answer" = "y" ]; then
s=$(echo "_$l")
fi
##-----------------------------------------------------------------------------

if [ ! -d results/22_CARD-AR_results"$s" ]; then
mkdir results/22_CARD-AR_results"$s"
fi

if [ ! -d results/22_CARD-AR_results"$s"/results ]; then
mkdir results/22_CARD-AR_results"$s"/results
fi

if [ ! -d results/22_CARD-AR_results"$s"/tmp ]; then
mkdir results/22_CARD-AR_results"$s"/tmp
fi

for F1 in $(cat $list); do

echo "runnning... CARD-ABR blast for $F1"

python2.7 /home/swapnil/pipeline/tools/rgi-3.2.0/rgi.py -t protein -i results/08_annotation/$F1/$F1.faa -o results/22_CARD-AR_results"$s"/results/$F1 > /dev/null 2>&1

rm results/22_CARD-AR_results"$s"/results/$F1.json

echo "finished... CARD-ABR blast for $F1"
done

#------------------------------------------------------------------------------
##calculate Abr-gene frequency

awk 'FNR>1' results/22_CARD-AR_results"$s"/results/*.txt > results/22_CARD-AR_results"$s"/tmp/all.txt

#-----

if grep -Fxq "ORF_ID	CONTIG	START	STOP	ORIENTATION	CUT_OFF	PASS_EVALUE	Best_Hit_evalue	Best_Hit_ARO	Best_Identities	ARO	ARO_name	Model_type	SNP	Best_Hit_ARO_category	ARO_category	PASS_bitscore	Best_Hit_bitscore	bit_score	Predicted_DNA	Predicted_Protein	CARD_Protein_Sequence	LABEL	ID	Model_id" results/22_CARD-AR_results"$s"/tmp/all.txt; then
:
else
ex -sc '1i|ORF_ID	CONTIG	START	STOP	ORIENTATION	CUT_OFF	PASS_EVALUE	Best_Hit_evalue	Best_Hit_ARO	Best_Identities	ARO	ARO_name	Model_type	SNP	Best_Hit_ARO_category	ARO_category	PASS_bitscore	Best_Hit_bitscore	bit_score	Predicted_DNA	Predicted_Protein	CARD_Protein_Sequence	LABEL	ID	Model_id' -cx results/22_CARD-AR_results"$s"/tmp/all.txt
fi

cp results/22_CARD-AR_results"$s"/tmp/all.txt results/22_CARD-AR_results"$s"/tmp/all.csv # create all.csv earlier and remove this step  
ssconvert results/22_CARD-AR_results"$s"/tmp/all.csv results/22_CARD-AR_results"$s"/tmp/all.csv.xlsx

#------
cp results/22_CARD-AR_results"$s"/tmp/all.txt results/22_CARD-AR_results"$s"/tmp/all.txt.1.tmp
sed -i 's/ /_/g' results/22_CARD-AR_results"$s"/tmp/all.txt.1.tmp
awk -F'\t' 'NR>1 {print $9,'\t',$15}' results/22_CARD-AR_results"$s"/tmp/all.txt.1.tmp > results/22_CARD-AR_results"$s"/tmp/all.txt.2.tmp
sed -i 's/  /===/' results/22_CARD-AR_results"$s"/tmp/all.txt.2.tmp 
cat results/22_CARD-AR_results"$s"/tmp/all.txt.2.tmp | sort | uniq > results/22_CARD-AR_results"$s"/tmp/all.txt.3.tmp 

rm results/22_CARD-AR_results"$s"/tmp/Abr_gene_frequency.csv.tmp
for F2 in $(cat results/22_CARD-AR_results"$s"/tmp/all.txt.3.tmp); do
V1=$(grep -c "$F2" results/22_CARD-AR_results"$s"/tmp/all.txt.2.tmp)
echo $V1 $F2 >> results/22_CARD-AR_results"$s"/tmp/Abr_gene_frequency.csv.tmp
done

cat results/22_CARD-AR_results"$s"/tmp/Abr_gene_frequency.csv.tmp | sort -k 1,1rn results/22_CARD-AR_results"$s"/tmp/Abr_gene_frequency.csv.tmp > results/22_CARD-AR_results"$s"/Abr_gene_frequency.csv

ex -sc '1i|frequency Abr-gene description' -cx results/22_CARD-AR_results"$s"/Abr_gene_frequency.csv

sed -i 's/===/ /g' results/22_CARD-AR_results"$s"/Abr_gene_frequency.csv
ssconvert results/22_CARD-AR_results"$s"/Abr_gene_frequency.csv results/22_CARD-AR_results"$s"/Abr_gene_frequency.csv.xlsx
rm results/22_CARD-AR_results"$s"/tmp/Abr_gene_frequency.csv.tmp

###############################################################################

echo "Completed.... CARD-ABR blast -------------------------------------------"

###############################################################################
## Creating 1-0 matrix of ABR results

echo "running matrix for ABR-DB"

for F1 in $(cat $list);do

awk -F'===' '{print $1}' results/22_CARD-AR_results"$s"/tmp/all.txt.3.tmp > results/22_CARD-AR_results"$s"/tmp/all.txt.4.tmp

awk -F'\t' '{print $9}' results/22_CARD-AR_results"$s"/results/$F1.txt | sed 's/ /_/g' > results/22_CARD-AR_results"$s"/tmp/$F1.gene-list-to-count.tmp1
	(rm results/22_CARD-AR_results"$s"/tmp/$F1.Abr-gene-count.tmp1) > /dev/null 2>&1 
	for V1 in $(cat results/22_CARD-AR_results"$s"/tmp/all.txt.4.tmp);do
	V2=$(awk '{count[$1]++} END {print count["'$V1'"]}' results/22_CARD-AR_results"$s"/tmp/$F1.gene-list-to-count.tmp1)
	echo $V2 >> results/22_CARD-AR_results"$s"/tmp/$F1.Abr-gene-count.tmp1
	done
awk '{for (i=1; i<= NF; i++) {if($i > 1) { $i=1; } } print }' results/22_CARD-AR_results"$s"/tmp/$F1.Abr-gene-count.tmp1 > results/22_CARD-AR_results"$s"/tmp/$F1.Abr-gene-count.tmp2
ex -sc '1i|'$F1'' -cx results/22_CARD-AR_results"$s"/tmp/$F1.Abr-gene-count.tmp1
ex -sc '1i|'$F1'' -cx results/22_CARD-AR_results"$s"/tmp/$F1.Abr-gene-count.tmp2
sed -i -e 's/^$/0/' results/22_CARD-AR_results"$s"/tmp/$F1.Abr-gene-count.tmp1
sed -i -e 's/^$/0/' results/22_CARD-AR_results"$s"/tmp/$F1.Abr-gene-count.tmp2
done

cp results/22_CARD-AR_results"$s"/tmp/all.txt.4.tmp results/22_CARD-AR_results"$s"/tmp/all.txt.4.2.tmp
sed -i 's/$/:c/' results/22_CARD-AR_results"$s"/tmp/all.txt.4.2.tmp
sed -i '1i\===\' results/22_CARD-AR_results"$s"/tmp/all.txt.4.2.tmp
paste results/22_CARD-AR_results"$s"/tmp/all.txt.4.2.tmp results/22_CARD-AR_results"$s"/tmp/*.Abr-gene-count.tmp1 > results/22_CARD-AR_results"$s"/tmp/matrix.csv.tmp1
paste results/22_CARD-AR_results"$s"/tmp/all.txt.4.2.tmp results/22_CARD-AR_results"$s"/tmp/*.Abr-gene-count.tmp2 > results/22_CARD-AR_results"$s"/tmp/matrix_1-0.csv.tmp2

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
}' results/22_CARD-AR_results"$s"/tmp/matrix.csv.tmp1 > results/22_CARD-AR_results"$s"/matrix.csv

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
}' results/22_CARD-AR_results"$s"/tmp/matrix_1-0.csv.tmp2 > results/22_CARD-AR_results"$s"/matrix_1-0.csv

sed -i 's/===//g' results/22_CARD-AR_results"$s"/matrix.csv
sed -i 's/===//g' results/22_CARD-AR_results"$s"/matrix_1-0.csv
sed -i 's/ /,/g' results/22_CARD-AR_results"$s"/matrix_1-0.csv 


rm results/21_blast_results/tmp/*.tmp1
rm results/21_blast_results/tmp/*.tmp2

echo "finished matrix for ABR-DB"

###############################################################################
