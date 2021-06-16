#!/bin/bash
###############################################################################
mkdir alignment_seperated_recent
mkdir tmp_recent
mkdir tmp_recent/fas
mkdir tmp_recent/fasta

##-----------------
perl split_multifasta.pl --input_file=core_gene_alignment.aln --output_dir=alignment_seperated_recent

awk 'NR>2 {print $4, $6}' recombinations_recent.txt | awk '!a[$1]++' > tmp_recent/lineage_strain.recent.tmp
awk 'NR>2 {print $4}' recombinations_recent.txt | awk '!a[$1]++' > tmp_recent/lineage.recent.tmp

##-----------------
for F1 in $(cat tmp_recent/lineage.recent.tmp); do
awk '$4 == "'$F1'" {print $1}' recombinations_recent.txt > tmp_recent/$F1.forward.lineage_information.recent.tmp
awk '$4 == "'$F1'" {print $2}' recombinations_recent.txt > tmp_recent/$F1.reverse.lineage_information.recent.tmp
done
##-----------------

for F1 in $(cat tmp_recent/lineage.recent.tmp); do
strain=$(awk '$1 == "'$F1'" {print $2}' tmp_recent/lineage_strain.recent.tmp)
samtools faidx alignment_seperated_recent/$strain.fsa

##-------
exec 3<tmp_recent/$F1.forward.lineage_information.recent.tmp
exec 4<tmp_recent/$F1.reverse.lineage_information.recent.tmp
while read forward <&3; read reverse <&4
do
samtools faidx alignment_seperated_recent/$strain.fsa "$strain":"$forward"-"$reverse" > tmp_recent/fas/$strain.$forward.$reverse.recent.fas
sed -i 's/-/N/g' tmp_recent/fas/$strain.$forward.$reverse.recent.fas
done
done
##----------------------------------------------------------
##----------------------------------------------------------

for F1 in $(cat tmp_recent/lineage.recent.tmp); do
strain=$(awk '$1 == "'$F1'" {print $2}' tmp_recent/lineage_strain.recent.tmp)
ls tmp_recent/fas/$strain.*.fas | sed 's/tmp_recent\/fas//g' | sed 's/\.fas//g' > tmp_recent/$strain.list.recent.tmp

cp /media/swapnil/Serratia/subproject_closed_chr/results/08_annotation/$strain/$strain.fsa tmp_recent/fasta/$strain.fsa
cp /media/swapnil/Serratia/subproject_closed_chr/results/08_annotation/$strain/$strain.gbk.csv tmp_recent/fasta/$strain.gbk.csv
sed -i 's/ /_/g' tmp_recent/fasta/$strain.gbk.csv
(makeblastdb -in tmp_recent/fasta/$strain.fsa -parse_seqids -dbtype nucl) > /dev/null 2>&1

##-------------
for F2 in $(cat tmp_recent/$strain.list.recent.tmp); do
V1=$(blastn -db tmp_recent/fasta/$strain.fsa -query tmp_recent/fas/$F2.fas -max_target_seqs 1 -perc_identity 100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk 'NR==1 {print $3, $4}')
if [ -z "$V1" ]; then
V1=$(blastn -db tmp_recent/fasta/$strain.fsa -query tmp_recent/fas/$F2.fas -max_target_seqs 1 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk 'NR==1 {print $3, $4}')
fi

if [ -z "$V1" ]; then
:
else
V2a=$(echo "$V1" | awk -F' ' '{print $1}')
V2=$(( $V2a - 200 ))
V3a=$(echo "$V1" | awk -F' ' '{print $2}')
V3=$(( $V3a + 200 ))
V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' tmp_recent/fasta/$strain.gbk.csv)
V5=$(echo "$V4" | awk -F' ' '{print $3}')
##----
	if [ -z "$V5" ]; then
	V2a=$(echo "$V1" | awk -F' ' '{print $1}')
	V2=$(( $V2a - 500 ))
	V3a=$(echo "$V1" | awk -F' ' '{print $2}')
	V3=$(( $V3a + 500 ))
	V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' tmp_recent/fasta/$strain.gbk.csv)
	V5=$(echo "$V4" | awk -F' ' '{print $3}')
	fi

	if [ -z "$V5" ]; then
	V2a=$(echo "$V1" | awk -F' ' '{print $1}')
	V2=$(( $V2a - 1000 ))
	V3a=$(echo "$V1" | awk -F' ' '{print $2}')
	V3=$(( $V3a + 1000 ))
	V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' tmp_recent/fasta/$strain.gbk.csv)
	V5=$(echo "$V4" | awk -F' ' '{print $3}')
	fi

	if [ -z "$V5" ]; then
	V2a=$(echo "$V1" | awk -F' ' '{print $1}')
	V2=$(( $V2a - 2000 ))
	V3a=$(echo "$V1" | awk -F' ' '{print $2}')
	V3=$(( $V3a + 2000 ))
	V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' tmp_recent/fasta/$strain.gbk.csv)
	V5=$(echo "$V4" | awk -F' ' '{print $3}')
	fi

	if [ -z "$V5" ]; then
	V2a=$(echo "$V1" | awk -F' ' '{print $1}')
	V2=$(( $V2a - 5000 ))
	V3a=$(echo "$V1" | awk -F' ' '{print $2}')
	V3=$(( $V3a + 5000 ))
	V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' tmp_recent/fasta/$strain.gbk.csv)
	V5=$(echo "$V4" | awk -F' ' '{print $3}')
	fi

	if [ -z "$V5" ]; then
	V4=$(echo "NA")
	fi
##----
S1=$(echo "$V4" | wc -l)
seq $S1 | xargs -I -- echo "$F2 $V2a $V3a $V2 $V3" >> tmp_recent/for_stat.1.csv
echo $V4 | awk -F' ' '{ for (i=3;i<=NF;i+=3) print $i }' | awk '{for (i=1;i<=NF;i++) print $i}' >> tmp_recent/for_stat.2.csv
echo $V4 | awk -F' ' '{ for (i=3;i<=NF;i+=3) print $i }' | awk '{for (i=1;i<=NF;i++) print $i}' >> tmp_recent/$strain.results.recent.tmp
fi
##----
done
paste tmp_recent/for_stat.1.csv tmp_recent/for_stat.2.csv > tmp_recent/for_stat.recent.csv
##-------------
## gene-accession and encoding proteins 

(cat tmp_recent/$strain.results.recent.tmp | sort -u ) > tmp_recent/$strain.results.gene-accession.recent.csv
for F3 in $(cat tmp_recent/$strain.results.gene-accession.recent.csv); do
V6=$(awk -F',' '{ if ($1 == "CDS" && $11 == "'$F3'") print $11, $15}' tmp_recent/fasta/$strain.gbk.csv)
V7=$(echo $V6 | awk '{print $2}')
if [ ${#V7} -le 2 ]; then 
V6=$(awk -F',' '{ if ($1 == "CDS" && $11 == "'$F3'") print $11, $16}' tmp_recent/fasta/$strain.gbk.csv)
	if [ ${#V6} -le 2 ]; then 
	V6=$(awk -F',' '{ if ($1 == "CDS" && $11 == "'$F3'") print $11, $17}' tmp_recent/fasta/$strain.gbk.csv)
	else
	:
	fi
echo $V6 >> $strain.results.recent.csv
else
echo $V6 >> $strain.results.recent.csv
fi

done
##-------------
done
##-------------
## adding lineages and bayesin factors to the gene-accession

for F1 in $(cat tmp_recent/lineage.recent.tmp); do
strain=$(awk '$1 == "'$F1'" {print $2}' tmp_recentl/lineage_strain.recent.tmp)

for F3 in $(cat tmp_recent/$strain.results.gene-accession.csv); do
T1=$(awk -F' ' '$6 == "'$F3'" {print $1}' tmp_recent/for_stat.recent.csv | awk -F'.' '{print $2, $3}' )
T2=$(echo $T1 | awk '{print $1}')
T3=$(echo $T1 | awk '{print $2}')
T4=$(awk '{ if ($1 == "'$T2'" && $2 == "'$T3'") print $3, $4, $5}' recombinations_recent.txt)
echo $F3 $T4 >> tmp_recent/$strain.for_stat.BF-added.recent.csv
paste tmp_recent/$strain.results.csv tmp_recent/$strain.for_stat.BF-added.recent.csv > $strain.results.final.recent.csv
done
sed -i '1s/^/gene_acc lng1 lng2 log(BF) gene_acc protein\n/' $strain.results.final.ancestral.csv
done
##-------------
done
##---------------------------


