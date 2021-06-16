#!/bin/bash
###############################################################################
#25 plasmid serach complete

###############################################################################

echo "started.... plasmid search complete -------------------------------------"

###############################################################################

if [ ! -d results/25_blast_plasmid_finder_complete_selected ]; then
mkdir results/25_blast_plasmid_finder_complete_selected
fi

#------------------------------------------------------------------------------
## Create database 
(rm data/reference-plasmids/all_ref_plasmids.fasta)> /dev/null 2>&1
cat data/reference-plasmids/*.fasta > data/reference-plasmids/all_ref_plasmids.fasta
(makeblastdb -in data/reference-plasmids/all_ref_plasmids.fasta -parse_seqids -dbtype nucl)> /dev/null 2>&1
echo "nucl blast database created for all_ref_plasmids.fasta"
#------------------------------------------------------------------------------

for F1 in $(cat list.txt)
do

if [ ! -d results/25_blast_plasmid_finder_complete_selected/$F1 ]; then
mkdir results/25_blast_plasmid_finder_complete_selected/$F1
fi

if [ ! -d results/25_blast_plasmid_finder_complete_selected/tmp ]; then
mkdir results/25_blast_plasmid_finder_complete_selected/tmp
fi

if [ ! -d results/25_blast_plasmid_finder_complete_selected/$F1/tmp ]; then
mkdir results/25_blast_plasmid_finder_complete_selected/$F1/tmp
fi

blastn -db data/reference-plasmids/all_ref_plasmids.fasta -query results/04_assembly/$F1/$F1.fasta -out results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.RefseqPlasmid.tmp -max_target_seqs 1 -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"

if grep -Fxq "sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.RefseqPlasmid.tmp; then
:
else
ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs' -cx results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.RefseqPlasmid.tmp
fi
sed -i '1s/ /\t/g' results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.RefseqPlasmid.tmp

echo "$F1" > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/strain_name.txt

#-------------------------------------------

grep ">" data/reference-plasmids/all_ref_plasmids.fasta | sed 's/>//g' | sort -k 1,1 -u > results/25_blast_plasmid_finder_complete_selected/tmp/plasmids.names.unique.tab

#-------------------------------------------

	for F2 in $(cat results/25_blast_plasmid_finder_complete_selected/tmp/plasmids.names.unique.tab);do
	echo "$F2" > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_accession_numbers_unique.tmp.tmp

	perl /home/swapnil/pipeline/tools/fastagrep.pl -f results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_accession_numbers_unique.tmp.tmp data/reference-plasmids/all_ref_plasmids.fasta > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F2.RefSeq.fasta

	awk '$1 == "'$F2'" {print $2}' results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.RefseqPlasmid.tmp > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.plasmid_accession_numbers_unique.nodes.tmp
	wc -l results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.plasmid_accession_numbers_unique.nodes.tmp | awk '{print $1}' > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.plasmid_accession_numbers_unique.nodes.count.tmp

	perl /home/swapnil/pipeline/tools/fastagrep.pl -f results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.plasmid_accession_numbers_unique.nodes.tmp results/04_assembly/$F1/$F1.fasta > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.fasta

#-------------------------------------------
## mapping contigs against reference

if [ -s results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.plasmid_accession_numbers_unique.nodes.tmp ]; then

cp results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F2.RefSeq.fasta /home/swapnil/pipeline/tools/mauve/
cp results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.fasta /home/swapnil/pipeline/tools/mauve/

WorDir=$(echo $PWD)

cd /home/swapnil/pipeline/tools/mauve
(java -Xmx5000m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output "$F1.$F2"_mauve -ref $F2.RefSeq.fasta  -draft $F1.$F2.fasta) > /dev/null 2>&1
rm $F2.RefSeq.fasta
rm $F1.$F2.fasta
cd $WorDir

mv /home/swapnil/pipeline/tools/mauve/"$F1.$F2"_mauve results/25_blast_plasmid_finder_complete_selected/$F1/tmp/
Var1=$(ls results/25_blast_plasmid_finder_complete_selected/$F1/tmp/"$F1.$F2"_mauve | sort -r | head -1)
cp results/25_blast_plasmid_finder_complete_selected/$F1/tmp/"$F1.$F2"_mauve/$Var1/$F1.$F2.fasta  results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.aligned.fasta
cp results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.aligned.fasta results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.aligned.joined.fasta
sed -i 's/>.*/NNNNNNNNNN/g' results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.aligned.joined.fasta
sed -i "1i "'>'$F1.$F2"" results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.aligned.joined.fasta
rm -r results/25_blast_plasmid_finder_complete_selected/$F1/tmp/"$F1.$F2"_mauve

# calculate ANI ---------------------------------------------------------------

mkdir /home/swapnil/pipeline/tools/output

(rm results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.tmp) > /dev/null 2>&1

(perl /home/swapnil/pipeline/tools/ANI_and_Total-aligned_plasmid_200.pl --fd /home/swapnil/pipeline/tools/blast-2.2.26/bin/formatdb --bl /home/swapnil/pipeline/tools/blast-2.2.26/bin/blastall --qr results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.aligned.joined.fasta --sb results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F2.RefSeq.fasta --od /home/swapnil/pipeline/tools/output >> results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.tmp) > /dev/null 2>&1

if [[ -s results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.tmp ]]; then
awk -F' ' '{print $2, $3, $4, $5}' results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.tmp > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.2.tmp
fi

awk -F' ' '{queryCov= 100*$2 / $3} END {print queryCov}' results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.2.tmp > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.query-covered.tmp ;

awk -F' ' '{refCov= 100*$2 / $4} END {print refCov}' results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.2.tmp > results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.ref-covered.tmp ;


rm -r /home/swapnil/pipeline/tools/output

#------------------------------------------------------------------------------
(rm results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_results.csv) > /dev/null 2>&1

paste results/25_blast_plasmid_finder_complete_selected/$F1/tmp/strain_name.txt results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_accession_numbers_unique.tmp.tmp results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.2.tmp results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.query-covered.tmp results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.ANI.ref-covered.tmp results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.plasmid_accession_numbers_unique.nodes.count.tmp >> results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_results.csv


# Annotation of plasmid -------------------------------------------------------

if [ ! -d results/08_annotation_plasmid ]; then
mkdir results/08_annotation_plasmid
fi

if [ ! -d results/08_annotation_plasmid/$F1 ]; then
mkdir results/08_annotation_plasmid/$F1
fi

/home/swapnil/tools/prokka-1.13/bin/prokka --quiet --outdir results/08_annotation_plasmid/$F1/$F1.$F2 --force --prefix $F1.$F2 --addgenes --locustag $F1.$F2 --strain $F1.$F2 --rnammer results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.$F2.aligned.joined.fasta

#------------------------------------------------------------------------------


cat results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_results.csv >> results/25_blast_plasmid_finder_complete_selected/all.csv
	
ex -sc '1i|strain plasmid ANI Aprox_alignment_length query_total_length subject_total_length %-Query-Covered %-Ref-covered number_of_contigs' -cx results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_results.csv

sed -i "s/ /\t/g" results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_results.csv
unoconv -i FilterOptions=09,,system,1 -f xls results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_results.csv
mv results/25_blast_plasmid_finder_complete_selected/$F1/tmp/$F1.plasmid_results.xls results/25_blast_plasmid_finder_complete_selected/$F1/

else
:
fi

echo "finished... Plasmid search for $F1"

done
done

ex -sc '1i|strain plasmid ANI Aprox_alignment_length query_total_length subject_total_length %-Query-Covered %-Ref-covered number_of_contigs' -cx results/25_blast_plasmid_finder_complete_selected/all.csv
sed -i "s/ /\t/g" results/25_blast_plasmid_finder_complete_selected/all.csv

unoconv -i FilterOptions=09,,system,1 -f xls results/25_blast_plasmid_finder_complete_selected/all.csv

###############################################################################
echo "Completed.. plasmid serach complete ------------------------------------"

###############################################################################
