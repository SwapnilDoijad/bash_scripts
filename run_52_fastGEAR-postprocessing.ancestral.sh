mkdir results/52_fastGEAR_ancestral
mkdir results/52_fastGEAR_ancestral/alignment_seperated_ancestral
mkdir results/52_fastGEAR_ancestral/tmp_ancestral
mkdir results/52_fastGEAR_ancestral/tmp_ancestral/fas
mkdir results/52_fastGEAR_ancestral/tmp_ancestral/fasta

##-----------------
perl /home/swapnil/pipeline/tools/split_multifasta.pl --input_file=results/52_fastGEAR/core_gene_alignment.aln --output_dir=results/52_fastGEAR_ancestral/alignment_seperated_ancestral

awk 'NR>1 {print $3, $4}' results/52_fastGEAR/output/lineage_information.txt | awk '!a[$1]++' > results/52_fastGEAR_ancestral/tmp_ancestral/lineage_strain.ancestral.tmp
awk 'NR>1 {print $3}' results/52_fastGEAR/output/lineage_information.txt | awk '!a[$1]++' > results/52_fastGEAR_ancestral/tmp_ancestral/lineage.ancestral.tmp

##-----------------
for F1 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/lineage.ancestral.tmp); do
awk '$4 == "'$F1'" {print $1}' results/52_fastGEAR/output/recombinations_ancestral.txt > results/52_fastGEAR_ancestral/tmp_ancestral/$F1.forward.lineage_information.ancestral.tmp
awk '$4 == "'$F1'" {print $2}' results/52_fastGEAR/output/recombinations_ancestral.txt > results/52_fastGEAR_ancestral/tmp_ancestral/$F1.reverse.lineage_information.ancestral.tmp
done
##---------------------------

for F1 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/lineage.ancestral.tmp); do
strain=$(awk '$1 == "'$F1'" {print $2}' results/52_fastGEAR_ancestral/tmp_ancestral/lineage_strain.ancestral.tmp)
samtools faidx results/52_fastGEAR_ancestral/alignment_seperated_ancestral/$strain.fsa

##-------
exec 3<results/52_fastGEAR_ancestral/tmp_ancestral/$F1.forward.lineage_information.ancestral.tmp
exec 4<results/52_fastGEAR_ancestral/tmp_ancestral/$F1.reverse.lineage_information.ancestral.tmp
while read forward <&3; read reverse <&4
do
samtools faidx results/52_fastGEAR_ancestral/alignment_seperated_ancestral/$strain.fsa "$strain":"$forward"-"$reverse" > results/52_fastGEAR_ancestral/tmp_ancestral/fas/$strain.$forward.$reverse.fas
sed -i 's/-/N/g' results/52_fastGEAR_ancestral/tmp_ancestral/fas/$strain.$forward.$reverse.fas
done
done
##----------------------------------------------------------
##----------------------------------------------------------

for F1 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/lineage.ancestral.tmp); do
strain=$(awk '$1 == "'$F1'" {print $2}' results/52_fastGEAR_ancestral/tmp_ancestral/lineage_strain.ancestral.tmp)
ls -1 results/52_fastGEAR_ancestral/tmp_ancestral/fas/$strain.*.fas | sed 's/results\/52_fastGEAR_ancestral\/tmp_ancestral\/fas\///g' | sed 's/\.fas//g' > results/52_fastGEAR_ancestral/tmp_ancestral/$strain.list.ancestral.tmp

cp /media/swapnil/Serratia/subproject_closed_chr/results/08_annotation/$strain/$strain.fsa results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.fsa
cp /media/swapnil/Serratia/subproject_closed_chr/results/08_annotation/$strain/$strain.gbk.csv results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv
sed -i 's/ /_/g' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv
(makeblastdb -in results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.fsa -parse_seqids -dbtype nucl) > /dev/null 2>&1

##-------------
for F2 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/$strain.list.ancestral.tmp); do
V1=$(blastn -db results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.fsa -query results/52_fastGEAR_ancestral/tmp_ancestral/fas/$F2.fas -max_target_seqs 1 -perc_identity 100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk 'NR==1 {print $3, $4}')
if [ -z "$V1" ]; then
V1=$(blastn -db results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.fsa -query results/52_fastGEAR_ancestral/tmp_ancestral/fas/$F2.fas -max_target_seqs 1 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk 'NR==1 {print $3, $4}')
fi

if [ -z "$V1" ]; then
:
else
V2a=$(echo "$V1" | awk -F' ' '{print $1}')
V2=$(( $V2a - 200 ))
V3a=$(echo "$V1" | awk -F' ' '{print $2}')
V3=$(( $V3a + 200 ))
V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv)
V5=$(echo "$V4" | awk -F' ' '{print $3}')
##----
	if [ -z "$V5" ]; then
	V2a=$(echo "$V1" | awk -F' ' '{print $1}')
	V2=$(( $V2a - 500 ))
	V3a=$(echo "$V1" | awk -F' ' '{print $2}')
	V3=$(( $V3a + 500 ))
	V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv)
	V5=$(echo "$V4" | awk -F' ' '{print $3}')
	fi

	if [ -z "$V5" ]; then
	V2a=$(echo "$V1" | awk -F' ' '{print $1}')
	V2=$(( $V2a - 1000 ))
	V3a=$(echo "$V1" | awk -F' ' '{print $2}')
	V3=$(( $V3a + 1000 ))
	V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv)
	V5=$(echo "$V4" | awk -F' ' '{print $3}')
	fi

	if [ -z "$V5" ]; then
	V2a=$(echo "$V1" | awk -F' ' '{print $1}')
	V2=$(( $V2a - 2000 ))
	V3a=$(echo "$V1" | awk -F' ' '{print $2}')
	V3=$(( $V3a + 2000 ))
	V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv)
	V5=$(echo "$V4" | awk -F' ' '{print $3}')
	fi

	if [ -z "$V5" ]; then
	V2a=$(echo "$V1" | awk -F' ' '{print $1}')
	V2=$(( $V2a - 5000 ))
	V3a=$(echo "$V1" | awk -F' ' '{print $2}')
	V3=$(( $V3a + 5000 ))
	V4=$(awk -F',' '{ if ($1 == "CDS" && $3 >= '$V2' && $4 <= '$V3') print $3, $4, $11}' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv)
	V5=$(echo "$V4" | awk -F' ' '{print $3}')
	fi

	if [ -z "$V5" ]; then
	V4=$(echo "NA")
	fi
##----
S1=$(echo "$V4" | wc -l)
seq $S1 | xargs -I -- echo "$F2 $V2a $V3a $V2 $V3" >> results/52_fastGEAR_ancestral/tmp_ancestral/for_stat.1.csv
echo $V4 | awk -F' ' '{ for (i=3;i<=NF;i+=3) print $i }' | awk '{for (i=1;i<=NF;i++) print $i}' >> results/52_fastGEAR_ancestral/tmp_ancestral/for_stat.2.csv
echo $V4 | awk -F' ' '{ for (i=3;i<=NF;i+=3) print $i }' | awk '{for (i=1;i<=NF;i++) print $i}' >> results/52_fastGEAR_ancestral/tmp_ancestral/$strain.results.ancestral.tmp
fi
##----
done
paste results/52_fastGEAR_ancestral/tmp_ancestral/for_stat.1.csv results/52_fastGEAR_ancestral/tmp_ancestral/for_stat.2.csv > results/52_fastGEAR_ancestral/tmp_ancestral/for_stat.ancestral.csv
##-------------
## gene-accession and encoding proteins 

(cat results/52_fastGEAR_ancestral/tmp_ancestral/$strain.results.ancestral.tmp | sort -u ) > results/52_fastGEAR_ancestral/tmp_ancestral/$strain.results.gene-accession.csv
for F3 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/$strain.results.gene-accession.csv); do
V6=$(awk -F',' '{ if ($1 == "CDS" && $11 == "'$F3'") print $11, $15}' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv)
V7=$(echo $V6 | awk '{print $2}')
if [ ${#V7} -le 2 ]; then 
V6=$(awk -F',' '{ if ($1 == "CDS" && $11 == "'$F3'") print $11, $16}' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv)
	if [ ${#V6} -le 2 ]; then 
	V6=$(awk -F',' '{ if ($1 == "CDS" && $11 == "'$F3'") print $11, $17}' results/52_fastGEAR_ancestral/tmp_ancestral/fasta/$strain.gbk.csv)
	else
	:
	fi
echo $V6 >> results/52_fastGEAR_ancestral/tmp_ancestral/$strain.results.csv
else
echo $V6 >> results/52_fastGEAR_ancestral/tmp_ancestral/$strain.results.csv
fi

done
##-------------
done
##-------------
## adding lineages and bayesin factors to the gene-accession

for F1 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/lineage.ancestral.tmp); do
strain=$(awk '$1 == "'$F1'" {print $2}' results/52_fastGEAR_ancestral/tmp_ancestral/lineage_strain.ancestral.tmp)

for F3 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/$strain.results.gene-accession.csv); do
T1=$(awk -F' ' '$6 == "'$F3'" {print $1}' results/52_fastGEAR_ancestral/tmp_ancestral/for_stat.ancestral.csv | awk -F'.' '{print $2, $3}' )
T2=$(echo $T1 | awk '{print $1}')
T3=$(echo $T1 | awk '{print $2}')
T4=$(awk '{ if ($1 == "'$T2'" && $2 == "'$T3'") print $3, $4, $5}' results/52_fastGEAR/output/recombinations_ancestral.txt)
echo $F3 $T4 >> results/52_fastGEAR_ancestral/tmp_ancestral/$strain.for_stat.BF-added.ancestral.csv
paste results/52_fastGEAR_ancestral/tmp_ancestral/$strain.for_stat.BF-added.ancestral.csv results/52_fastGEAR_ancestral/tmp_ancestral/$strain.results.csv > results/52_fastGEAR_ancestral/$strain.results.final.ancestral.1.csv
done
sed '1s/^/gene_acc lng1 lng2 log(BF) gene_acc protein\n/' results/52_fastGEAR_ancestral/$strain.results.final.ancestral.1.csv > results/52_fastGEAR_ancestral/$strain.results.final.ancestral.csv
done
##-------------
done
##---------------------------

cat results/52_fastGEAR_ancestral/*.results.final.ancestral.1.csv >> results/52_fastGEAR_ancestral/all.results.final.ancestral.csv
awk '{print $1}' results/52_fastGEAR_ancestral/all.results.final.ancestral.csv | sed 's/......$//' | sort -u > results/52_fastGEAR_ancestral/tmp_ancestral/strain.post-processing.ancestral.tmp

for F4 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/strain.post-processing.ancestral.tmp); do
awk 'NR>1 {print $1}' results/52_fastGEAR_ancestral/$F4.results.final.ancestral.csv > results/52_fastGEAR_ancestral/tmp/$F4.recombined_genes.csv
done

for F5 in $(cat results/52_fastGEAR_ancestral/tmp_ancestral/strain.post-processing.ancestral.tmp); do
perl  /home/swapnil/pipeline/tools/fastagrep.pl -f results/52_fastGEAR_ancestral/tmp/$F5.recombined_genes.csv results/08_annotation/raw_files/$F5/$F5.faa > /media/swapnil/Serratia/subproject_closed_chr/results/52_fastGEAR_ancestral/tmp/$F5.recombined_genes.csv.faa
done

sed -i 's/ .*//g' results/52_fastGEAR_ancestral/tmp/*.recombined_genes.csv.faa
cat results/52_fastGEAR_ancestral/tmp/*.recombined_genes.csv.faa > results/52_fastGEAR_ancestral/tmp/all.recombined_genes.csv.faa

