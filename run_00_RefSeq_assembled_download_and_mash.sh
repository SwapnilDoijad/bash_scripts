#!/bin/bash
###############################################################################
# 00 harvest suite

###############################################################################

echo "RefSeq genome downloading ----------------------------------------------"

###############################################################################
if [ ! -d results ]; then
mkdir results
fi

if [ ! -d results/00_ref ]; then
mkdir results/00_ref
fi

if [ ! -d results/00_ref/RefSeq ]; then
mkdir results/00_ref/RefSeq
fi

if [ ! -d results/00_ref/RefSeq/RefSeq_fasta ]; then
mkdir results/00_ref/RefSeq/RefSeq_fasta
fi

if [ ! -d results/00_ref/RefSeq/RefSeq_msh ]; then
mkdir results/00_ref/RefSeq/RefSeq_msh
fi

if [ ! -d results/00_ref/RefSeq/tmp ]; then
mkdir results/00_ref/RefSeq/tmp
fi
 
awk -F'\t' '{print $19}' genomes_proks.txt | sed '/ftp/!d' > results/00_ref/RefSeq/tmp/refseq.list.1.tmp
awk -F'/' '{print $10}' results/00_ref/RefSeq/tmp/refseq.list.1.tmp > results/00_ref/RefSeq/tmp/refseq.list.2.tmp
sed -i 's/$/\//g' results/00_ref/RefSeq/tmp/refseq.list.1.tmp
paste results/00_ref/RefSeq/tmp/refseq.list.1.tmp results/00_ref/RefSeq/tmp/refseq.list.2.tmp > results/00_ref/RefSeq/tmp/refseq.download-link.list.tmp
sed -i 's/\t//g' results/00_ref/RefSeq/tmp/refseq.download-link.list.tmp

for F1 in $(cat results/00_ref/RefSeq/tmp/refseq.download-link.list.tmp);do
	V1=$(echo $F1 | awk -F'/' '{print $10}')
	echo "Downloading $V1"
(wget -P results/00_ref/RefSeq/tmp "$F1"_genomic.fna.gz) > /dev/null 2>&1 ;

gunzip results/00_ref/RefSeq/tmp/"$V1"_genomic.fna.gz

V2=$(head -n 1 results/00_ref/RefSeq/tmp/"$V1"_genomic.fna | sed 's/\..*$//' | sed 's/>//g')



mv results/00_ref/RefSeq/tmp/"$V1"_genomic.fna results/00_ref/RefSeq/tmp/"$V2".fasta 

#----------------------------------

V3=$(grep -c ">" results/00_ref/RefSeq/tmp/"$V2".fasta)
echo $V3 > results/00_ref/RefSeq/tmp/$V2.all-contigs-number.tmp
if [ "$V3" -lt 1000 ];then
	cp results/00_ref/RefSeq/tmp/"$V2".fasta results/00_ref/RefSeq/RefSeq_fasta/

	echo $V2 >> results/00_ref/RefSeq/RefSeq.final.list.txt

#----------------------------------
	## Mash distance
	cp results/00_ref/RefSeq/RefSeq_fasta/$V2.fasta  results/00_ref/RefSeq/tmp/$V2.fasta.tmp
	sed -i 's/>/>gi|00000000|XXX|00000000| /g' results/00_ref/RefSeq/tmp/$V2.fasta.tmp
	echo sketching $V2 
	(mash sketch -o results/00_ref/RefSeq/RefSeq_msh/$V2.msh results/00_ref/RefSeq/tmp/$V2.fasta.tmp)> /dev/null 2>&1
	rm results/00_ref/RefSeq/tmp/$V2.fasta.tmp
else
:
fi

done

###############################################################################

echo "RefSeq genome downloading completed ------------------------------------"

###############################################################################




