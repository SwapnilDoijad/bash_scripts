#!/bin/bash
###############################################################################
#24 plasmid finder
###############################################################################
echo "started.... Plasmid search ---------------------------------------------"
###############################################################################
echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")
#------------------------------------------------------------------------------

if [ ! -d results/24_blast_plasmid_finder_"$l" ]; then
mkdir results/24_blast_plasmid_finder_"$l"
fi

#------------------------------------------------------------------------------

for F1 in $(cat list.txt)
do

echo "runnning... Plasmid search for $F1"

mkdir results/24_blast_plasmid_finder_"$l"/$F1

blastn -db /home/swapnil/pipeline/tools/databases/plasmid/plasmids.fasta -query results/04_assembly/raw_files/$F1/$F1.fasta -out results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp -max_target_seqs 1 -evalue 1e-100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
#------------------------------------------------------------------------------

if [[ -s results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp ]] ; then

#------------------------------------------------------------------------------

cp results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp2
sed -i 's/NODE_/'$F1'_NODE_/g' results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp2
awk '{print $2}' results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp > results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp3

#------------------------------------------------------------------------------
#isolate plasmid fasta and write to folder

/home/swapnil/pipeline/tools/fastagrep.pl -f results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp3 results/04_assembly/raw_files/$F1/$F1.fasta > results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_contig.all.fasta

split_fasta.pl -prefix $F1. results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_contig.all.fasta


while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=$F1.plasmid_contig.${line#>}.fa
        echo $line > results/24_blast_plasmid_finder_"$l"/$F1/$outfile
	echo $outfile >> results/24_blast_plasmid_finder_"$l"/$F1/plasmid-contigs.list ; chmod a+x results/24_blast_plasmid_finder_"$l"/$F1/plasmid-contigs.list
    else
        echo $line >> results/24_blast_plasmid_finder_"$l"/$F1/$outfile
    fi
done < results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_contig.all.fasta
#------------------------------------------------------------------------------

cat results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp2 >> results/24_blast_plasmid_finder_"$l"/plasmid_finder_blast.tab


#------------------------------------------------------------------------------
# copy contigs adn label it with plasmid-replicon

if [ -f results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp3 ]; then
cp results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp.tmp
sed -i 's/\t/===/g' results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp.tmp

for V1 in $(cat  results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_finder.tmp.tmp); do
V2=$(echo $V1 | awk -F'_' '{print $1}')
V3=$(echo $V1 | awk -F'===' '{print $2}')
cp results/24_blast_plasmid_finder_"$l"/$F1/$F1.plasmid_contig.$V3.fa results/24_blast_plasmid_finder_"$l"/$F1/$F1.$V2.$V3.fasta
done

fi

#------------------------------------------------------------------------------
#BLAST contigs identified as of plasmids against refSeq-Plasmid database

for F2 in $(cat results/24_blast_plasmid_finder_"$l"/$F1/plasmid-contigs.list)
do
blastn -db /home/swapnil/pipeline/tools/databases/plasmid/refseq_plasmid.fasta -query results/24_blast_plasmid_finder_"$l"/$F1/$F2 -out results/24_blast_plasmid_finder_"$l"/$F1/refseq.$F2.tmp -max_target_seqs 1 -evalue 1e-100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
sed -i 's/NODE_/'$F1'_NODE_/g' results/24_blast_plasmid_finder_"$l"/$F1/refseq.$F2.tmp
head -n 1 results/24_blast_plasmid_finder_"$l"/$F1/refseq.$F2.tmp >> results/24_blast_plasmid_finder_"$l"/$F1/refseq.all.tmp
done

#------------------------------------------------------------------------------

awk -f /home/swapnil/pipeline/tools/vlookup-plasmid-refseq.awk /home/swapnil/pipeline/tools/databases/plasmid/refseq_plasmid.fasta.annotation.tab results/24_blast_plasmid_finder_"$l"/$F1/refseq.all.tmp > results/24_blast_plasmid_finder_"$l"/$F1/vlookup.tmp
paste results/24_blast_plasmid_finder_"$l"/$F1/refseq.all.tmp results/24_blast_plasmid_finder_"$l"/$F1/vlookup.tmp >> results/24_blast_plasmid_finder_"$l"/RefSeq-plasmids-blast.csv

#------------------------------------------------------------------------------
else
:
fi
#------------------------------------------------------------------------------
echo "finished... Plasmid search for $F1"
done

if grep -Fxq "sseqid	qseqid	sstart	send	qstart	qend	slen	qlen	evalue	bitscore	length	mismatch	gaps	pident	qcovs" results/24_blast_plasmid_finder_"$l"/plasmid_finder_blast.tab; then
:
else
ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs' -cx results/24_blast_plasmid_finder_"$l"/plasmid_finder_blast.tab
fi
sed -i '1s/ /\t/g' results/24_blast_plasmid_finder_"$l"/plasmid_finder_blast.tab


if grep -Fxq "sseqid	qseqid	sstart	send	qstart	qend	slen	qlen	evalue	bitscore	length	mismatch	gaps	pident	qcovs	plasmid" results/24_blast_plasmid_finder_"$l"/RefSeq-plasmids-blast.csv; then
:
else
ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs plasmid' -cx results/24_blast_plasmid_finder_"$l"/RefSeq-plasmids-blast.csv
fi
sed -i 's/===/ /g' results/24_blast_plasmid_finder_"$l"/RefSeq-plasmids-blast.csv
###############################################################################
# copy plasmid fasta as per replicon

if [ ! -d results/24_blast_plasmid_finder_"$l"/all_fasta_as_per_replicon ]; then
mkdir results/24_blast_plasmid_finder_"$l"/all_fasta_as_per_replicon
fi

awk -F'_' 'FNR >1 {print $1}' results/24_blast_plasmid_finder_"$l"/plasmid_finder_blast.tab | sort -u > results/24_blast_plasmid_finder_"$l"/tmp/sorted_plasmid_all_fasta_as_per_replicon.txt

for F3 in $(cat results/24_blast_plasmid_finder_"$l"/tmp/sorted_plasmid_all_fasta_as_per_replicon.txt);do

if [ ! -d results/24_blast_plasmid_finder_"$l"/all_fasta_as_per_replicon/$F3 ]; then
mkdir results/24_blast_plasmid_finder_"$l"/all_fasta_as_per_replicon/$F3
fi

cp results/24_blast_plasmid_finder_"$l"/results/*/*.$F3.*.fasta results/24_blast_plasmid_finder_"$l"/all_fasta_as_per_replicon/$F3/

done

###############################################################################

echo "Completed.. Plasmid search ---------------------------------------------"

###############################################################################
