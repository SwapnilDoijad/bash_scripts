#!/bin/bash
###############################################################################
echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

echo "want to suffix files? type "y" else press enter to continue"
read answer
if [ "$answer" = "y" ]; then
s=$(echo "_$l")
fi

(mkdir results/17_roary"$s"/truncation) > /dev/null 2>&1
(mkdir results/17_roary"$s"/truncation/tmp) > /dev/null 2>&1

awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", "_", $i) } 1' results/17_roary"$s"/roary_results/gene_presence_absence.csv > results/17_roary"$s"/truncation/gene_presence_absence.1.csv ;

awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(" ", "_", $i) } 1' results/17_roary"$s"/truncation/gene_presence_absence.1.csv > results/17_roary"$s"/truncation/gene_presence_absence.csv ;

rm results/17_roary"$s"/truncation/gene_presence_absence.1.csv ;

head -1 results/17_roary"$s"/truncation/gene_presence_absence.csv | sed "s/\r//" | awk -F',' '{out=""; for(i=15;i<=NF;i++){out=out"\n"$i}; print out}' > results/17_roary"$s"/truncation/tmp/isolate.list

awk -F',' '($14 - $12) > 30  {print $0}' results/17_roary"$s"/truncation/gene_presence_absence.csv > results/17_roary"$s"/truncation/tmp/difference.tmp #loss of 10 A.A. 

awk -F',' '{print $1}' results/17_roary"$s"/truncation/tmp/difference.tmp | sed 's/-/_/g' | sed "s/'/_/g" > results/17_roary"$s"/truncation/tmp/difference.list.tmp

##-----------------------------------------------------------------------------
(rm results/17_roary"$s"/truncation/truncated_gene_list_and_true_length.csv) > /dev/null 2>&1
for trun_gene in $(cat results/17_roary"$s"/truncation/tmp/difference.list.tmp); do

sed 's/-//g' results/17_roary"$s"/roary_results/pan_genome_sequences/$trun_gene.fa.aln | bioawk -c fastx '{ print $name, length($seq) }' | sort -n -k2 > results/17_roary"$s"/truncation/tmp/$trun_gene.sorted.tmp

##------------------------------------- 
## find gene true length
gene_true_length=$(sed 's/-//g' results/17_roary"$s"/roary_results/pan_genome_sequences/$trun_gene.fa.aln | bioawk -c fastx '{ print $name, length($seq) }' | sort -n -k2 | awk '{print $2}' | uniq -c | sort -n -r -k2 | awk 'FNR==1 {print $2}')
if [ "$gene_true_length" -gt 3 ]; then
gene_true_length=$(sed 's/-//g' results/17_roary"$s"/roary_results/pan_genome_sequences/$trun_gene.fa.aln | bioawk -c fastx '{ print $name, length($seq) }' | sort -n -k2 | awk '{print $2}' | uniq -c | sort -n -r -k2 | awk 'FNR==1 {print $2}')
else
gene_true_length=$(sed 's/-//g' results/17_roary"$s"/roary_results/pan_genome_sequences/$trun_gene.fa.aln | bioawk -c fastx '{ print $name, length($seq) }' | sort -n -k2 | awk '{print $2}' | uniq -c | sort -n -r -k2 | awk 'FNR==2 {print $2}')
	if [ "$gene_true_length" -gt 3 ]; then
	gene_true_length=$(sed 's/-//g' results/17_roary"$s"/roary_results/pan_genome_sequences/$trun_gene.fa.aln | bioawk -c fastx '{ print $name, length($seq) }' | sort -n -k2 | awk '{print $2}' | uniq -c | sort -n -r -k2 | awk 'FNR==1 {print $2}')
	else
	gene_true_length=$(sed 's/-//g' results/17_roary"$s"/roary_results/pan_genome_sequences/$trun_gene.fa.aln | bioawk -c fastx '{ print $name, length($seq) }' | sort -n -k2 | awk '{print $2}' | uniq -c | sort -n -r -k2 | awk 'FNR==3 {print $2}')
	fi
fi
##------------------------------------- 

echo $trun_gene $gene_true_length >> results/17_roary"$s"/truncation/truncated_gene_list_and_true_length.csv

awk -F'\t' '$2 < '$gene_true_length' {print $1, $2}' results/17_roary"$s"/truncation/tmp/$trun_gene.sorted.tmp > results/17_roary"$s"/truncation/tmp/$trun_gene.truncated.list

#----------------------
(rm results/17_roary"$s"/truncation/tmp/$trun_gene.matrix.tmp) > /dev/null 2>&1
echo "$trun_gene" > results/17_roary"$s"/truncation/tmp/$trun_gene.matrix.tmp
for isolate in $(cat results/17_roary"$s"/truncation/tmp/isolate.list); do
V2=$(grep $isolate results/17_roary"$s"/truncation/tmp/$trun_gene.truncated.list)
if [ -z "$V2" ]; then
echo "1" >> results/17_roary"$s"/truncation/tmp/$trun_gene.matrix.tmp
else
echo "0" >> results/17_roary"$s"/truncation/tmp/$trun_gene.matrix.tmp
fi
done
##----------------------
# gen
(rm results/17_roary"$s"/truncation/tmp/$trun_gene.trunc_mtrx.tmp) > /dev/null 2>&1
echo "$trun_gene" > results/17_roary"$s"/truncation/tmp/$trun_gene.trunc_mtrx.tmp
for isolate in $(cat results/17_roary"$s"/truncation/tmp/isolate.list); do
gene_length=$( grep $isolate results/17_roary"$s"/truncation/tmp/$trun_gene.sorted.tmp | awk '{print $2}' | sort -r -n | head -1);
if [ -z "$gene_length" ]; then
echo "$isolate $trun_gene absent" >> test.trunc_mtrx.tmp
echo "0" >> results/17_roary"$s"/truncation/tmp/$trun_gene.trunc_mtrx.tmp
else
V3=$(echo "($gene_length*100)/$gene_true_length" | bc )
echo "$isolate $trun_gene $V3" >> test.trunc_mtrx.tmp
echo "$V3" >> results/17_roary"$s"/truncation/tmp/$trun_gene.trunc_mtrx.tmp
fi
done
#----------------------

done

##-----------------------------------------------------------------------------
## creating matrix
cp results/17_roary"$s"/truncation/tmp/isolate.list results/17_roary"$s"/truncation/tmp/00.isolate.matrix.tmp

#paste results/17_roary"$s"/truncation/tmp/00.isolate.matrix.tmp results/17_roary"$s"/truncation/tmp/*.matrix.tmp > results/17_roary"$s"/truncation/matrix.csv

# distance-matix-format-for_1-0_format
ls -1 results/17_roary"$s"/truncation/tmp/*.matrix.tmp  | split -l 500 -d - lists
for list in lists*; do
paste $(cat $list) > results/17_roary"$s"/truncation/merge.${list##lists}; 
done
paste results/17_roary"$s"/truncation/merge.* > results/17_roary"$s"/truncation/tmp/all.1-0.csv
rm results/17_roary"$s"/truncation/merge.* 
rm lists*


# distance-matix-format-for_truncation_format
cp results/17_roary"$s"/truncation/tmp/isolate.list results/17_roary"$s"/truncation/tmp/00.trunc_mtrx.tmp ;

ls -1 results/17_roary"$s"/truncation/tmp/*.trunc_mtrx.tmp  | split -l 500 -d - lists
for list in lists*; do
paste $(cat $list) > results/17_roary"$s"/truncation/merge.${list##lists}; 
done
paste results/17_roary"$s"/truncation/merge.* > results/17_roary"$s"/truncation/tmp/all.trunc.csv
rm results/17_roary"$s"/truncation/merge.* 
rm lists*

##-----------------------------------------------------------------------------
bash /home/swapnil/pipeline/tools/transpose.roary_truncation.sh 
