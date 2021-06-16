#!/bin/bash
###############################################################################
#32 variscan

#codeml for dN/dS not workiong properly, also taking time

###############################################################################
echo "started.... variscan ---------------------------------------------------"
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

if [ ! -d results/32_variscan"$s" ]; then
mkdir results/32_variscan"$s"
fi

if [ ! -d results/32_variscan"$s"/tmp ]; then
mkdir results/32_variscan"$s"/tmp
fi

if [ ! -d results/32_variscan"$s"/codeml ]; then
mkdir results/32_variscan"$s"/codeml
fi

if [ ! -d results/32_variscan"$s"/codeml/tmp ]; then
mkdir results/32_variscan"$s"/codeml/tmp
fi

cp /home/swapnil/pipeline/tools/codeml.ctl  /home/swapnil/

ls results/17_roary"$s"/roary_results/pan_genome_sequences/*.fa.aln | head -2 | awk -F'/' '{print $5}' | awk -F'.' '{print $1}' > results/32_variscan"$s"/tmp/list-gene-variscan.txt
###############################################################################

(rm results/32_variscan"$s"/tmp/all.results.txt) > /dev/null 2>&1

for gene in $(cat results/32_variscan"$s"/tmp/list-gene-variscan.txt);do
echo "running $gene"
##----------------------------------------------------------------------------
## calcualting dN/dS aka Ka/Ks ratio & Ts/Tv ratio

#(raxmlHPC -m GTRGAMMA -p 12345 -s results/17_roary"$s"/roary_results/pan_genome_sequences/$gene.fa.aln -n raxml) > /dev/null 2>&1
#(rm results/17_roary"$s"/roary_results/pan_genome_sequences/$gene.fa.aln.reduced) > /dev/null 2>&1
#(cp RAxML_bestTree.raxml /home/swapnil/tmp.tree) > /dev/null 2>&1
#(cp results/17_roary"$s"/roary_results/pan_genome_sequences/$gene.fa.aln /home/swapnil/tmp.aln) > /dev/null 2>&1
#(rm *.raxml) > /dev/null 2>&1

#Cur_dir=$(pwd)
#cd  /home/swapnil/
#(
#codeml
#) > /dev/null 2>&1 ;
#kappa=$(grep "kappa" tmp | awk '{print $4}')
#omega=$(grep "omega" tmp | awk '{print $4}')
#if [ -z "$kappa" ] ; then
#kappa="0"
#fi

#if [ -z "$omega" ] ; then
#omega="0"
#fi

#echo $gene $kappa $omega > $Cur_dir/results/32_variscan"$s"/codeml/tmp/$gene.txt
#rm 2NG.dN
#rm 2NG.dS
#rm 2NG.t
#rm 4fold.nuc
#rm lnf
#rm rst
#rm rst1
#rm rub
#rm tmp
#cd $Cur_dir

##----------------------------------------------------------------------------
cp results/17_roary"$s"/roary_results/pan_genome_sequences/$gene.fa.aln results/tmp.fa.aln

python2.7 /home/swapnil/pipeline/tools/convert_sequence_format.py

mv results/tmp.phy results/32_variscan"$s"/tmp/$gene.phy

(variscan results/32_variscan"$s"/tmp/$gene.phy /home/swapnil/pipeline/tools/variscan_conf.conf > results/32_variscan"$s"/tmp/$gene.results.txt) > /dev/null 2>&1

# Number-of-sites-analyzed-in-gene
V1=$(grep '#Num of analysed sites: ' results/32_variscan"$s"/tmp/$gene.results.txt | sed 's/^.*: //')
if [ -z "$V1" ] ; then
V1a="0"
else
V1a=$V1
fi

# gene present-in
V2=$(grep -c ">" results/17_roary"$s"/roary_results/pan_genome_sequences/$gene.fa.aln)
if [ -z "$V2" ] ; then
V2a="0"
else
V2a=$V2
fi

# Number-of-Polymorphic-sites(S)
V4=$(tail -1 results/32_variscan"$s"/tmp/$gene.results.txt | awk '{print $2}')
if [ -z "$V4" ] ; then
V4a="0"
else
V4a=$V4
fi

# Total-number-of-mutations(Eta)
V5=$(tail -1 results/32_variscan"$s"/tmp/$gene.results.txt | awk '{print $3}')
if [ -z "$V5" ] ; then
V5a="0"
else
V5a=$V5
fi

(printf "%.2f" $(($V4*100/$V1)) | sed 's/,/\./g' > results/32_variscan"$s"/tmp/$gene.1.tmp) > /dev/null 2>&1
(printf "%.2f" $(($V5*100/$V1)) | sed 's/,/\./g' > results/32_variscan"$s"/tmp/$gene.2.tmp) > /dev/null 2>&1

# S% 
V6=$(cat results/32_variscan"$s"/tmp/$gene.1.tmp)
if [ -z "$V6" ] ; then
V6a="0"
else
V6a=$V6
fi

# Eta%
V7=$(cat results/32_variscan"$s"/tmp/$gene.2.tmp)
if [ -z "$V7" ] ; then
V7a="0"
else
V7a=$V7
fi

(awk "BEGIN {print $V5a/$V4a}" > results/32_variscan"$s"/tmp/$gene.3.tmp) > /dev/null 2>&1
if [ -f results/32_variscan"$s"/tmp/$gene.3.tmp ]; then
V8=$(cat results/32_variscan"$s"/tmp/$gene.3.tmp)
	if [ -z "$V8" ] ; then
	V8a="0"
	else
	V8a=$V8
	fi
else
V8a="0"
fi

# Pi
V9=$(tail -1 results/32_variscan"$s"/tmp/$gene.results.txt | awk '{print $5}')
if [ -z "$V9" ] ; then
V9a="0"
else
V9a=$V9
fi

# Theta
V10=$(tail -1 results/32_variscan"$s"/tmp/$gene.results.txt | awk '{print $6}')
if [ -z "$V10" ] ; then
V10a="0"
else
V10a=$V10
fi

# Tajima_D
V11=$(tail -1 results/32_variscan"$s"/tmp/$gene.results.txt | awk '{print $7}')
if [ -z "$V11" ] ; then
V11a="0"
else
V11a=$V11
fi

# FuLi_Dstar
V12=$(tail -1 results/32_variscan"$s"/tmp/$gene.results.txt | awk '{print $8}')
if [ -z "$V12" ] ; then
V12a="0"
else
V12a=$V12
fi

# FuLi_Fstar
V13=$(tail -1 results/32_variscan"$s"/tmp/$gene.results.txt | awk '{print $9}')
if [ -z "$V13" ] ; then
V13a="0"
else
V13a=$V13
fi

rm results/32_variscan"$s"/tmp/$gene.1.tmp
rm results/32_variscan"$s"/tmp/$gene.2.tmp
rm results/32_variscan"$s"/tmp/$gene.3.tmp

echo $gene $kappa $omega $V2a $V1a $V4a $V6a $V5a $V7a $V8a $V9a $V10a $V11a $V12a $V13a >> results/32_variscan"$s"/tmp/all.results.txt

#rm results/32_variscan"$s"/tmp/$gene.phy
#rm results/32_variscan"$s"/tmp/$gene.results.txt

done

rm results/32_variscan"$s"/tmp/list-gene-variscan.txt

cp results/17_roary"$s"/roary_results/gene_presence_absence.csv results/32_variscan"$s"/tmp/

sed -i 's/ /_/g' results/32_variscan"$s"/tmp/gene_presence_absence.csv
sed -i 's/,/_/g' results/32_variscan"$s"/tmp/gene_presence_absence.csv
sed -i 's/"_"/===/g' results/32_variscan"$s"/tmp/gene_presence_absence.csv

awk -F"===" '{print $1, $3}' results/32_variscan"$s"/tmp/gene_presence_absence.csv > results/32_variscan"$s"/tmp/tmp.tmp

sed -i 's/\"//g' results/32_variscan"$s"/tmp/tmp.tmp

awk -f /home/swapnil/pipeline/tools/vlookup_variscan.awk results/32_variscan"$s"/tmp/tmp.tmp results/32_variscan"$s"/tmp/all.results.txt > results/32_variscan"$s"/tmp/tmp2.tmp

paste results/32_variscan"$s"/tmp/all.results.txt results/32_variscan"$s"/tmp/tmp2.tmp > results/32_variscan"$s"/tmp/results.tmp

sort -t' ' -k5,5rn -k7,7rn results/32_variscan"$s"/tmp/results.tmp > results/32_variscan"$s"/tmp/results.csv

sed -i '1 i\gene kappa omega present-in Number-of-sites-analyzed-in-gene Number-of-Polymorphic-sites(S) S% Total-number-of-mutations(Eta) Eta% Eta/S-ratio Pi Theta Tajima_D FuLi_Dstar FuLi_Fstar gene-product' results/32_variscan"$s"/tmp/results.csv

sed -i 's/ /\t/g' results/32_variscan"$s"/tmp/results.csv

ssconvert results/32_variscan"$s"/tmp/results.csv results/32_variscan"$s"/results.xlsx

rm results/tmp.fa.aln
rm /home/swapnil/tmp.tree
rm /home/swapnil/tmp.aln
rm /home/swapnil/codeml.ctl
###############################################################################

echo "Completed.. variscan ---------------------------------------------------"

###############################################################################


