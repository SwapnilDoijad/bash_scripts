#!/bin/bash
###############################################################################
#25 plasmid serach complete

# makeblastdb -in 310519_all.fsa -parse_seqids -dbtype nucl

###############################################################################

echo "started.... plasmid search complete -------------------------------------"

###############################################################################

echo "provide list file (for e.g. all)"
read l
s=$(echo "_$l")
list=$(echo "list.$l.txt")
echo "provide plasmid accession number (for e.g. NC_021845.1)"
read sp_p

#------------------------------------------------------------------------------
(mkdir results/25_blast_plasmid_finder_complete"$s")> /dev/null 2>&1
(mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/)> /dev/null 2>&1
(mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/)> /dev/null 2>&1


#------------------------------------------------------------------------------
echo "Isolate Plasmid contigs_involved size(bp) Ref_size(bp) %covered pCov gCov replicon WGS_replicons contigs" > results/25_blast_plasmid_finder_complete"$s"/$sp_p/results.csv

(rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$WGS_replicon.list)> /dev/null 2>&1
for F1 in $(cat $list); do

(mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1)> /dev/null 2>&1
(mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp)> /dev/null 2>&1
(mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/WGS)> /dev/null 2>&1
(mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta)> /dev/null 2>&1
(mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/plasmid_replicon)> /dev/null 2>&1

echo "runnning... Plasmid search for $F1"

##---------------------------
## coverage calculations

grep ">" results/04_assembly/raw_files/$F1/$F1.fasta > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.strain_genome_coverage-calculation.tmp

genome_cov=$(awk -F '_' '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.strain_genome_coverage-calculation.tmp)

##---------------------------
echo "getting WGS-replicon for $F1"
blastn -db /home/swapnil/pipeline/tools/databases/plasmid_replicons/310519_all.fsa -query results/04_assembly/raw_files/$F1/$F1.fasta -out results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/WGS/$F1.blast.tmp -max_target_seqs 1 -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
WGS_replicon=$(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/WGS/$F1.blast.tmp | awk '{print $1}' | sed '/^$/d' | tr "\n" "_" )
if [ -z "$WGS_replicon" ]; then 
WGS_replicon=$(echo "NA")
else
cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/WGS/$F1.blast.tmp | awk '{print $2}' | sed -e "s/$/\.fa/" > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.WGS_replicon.contigs.list
fi
echo $F1 $WGS_replicon > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.WGS_replicon.list

##---------------------------
## fastgrep > 5000 bp length contigs and split assembly to seperate contig files
cp results/04_assembly/raw_files/$F1/$F1.fasta results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp
grep -F ">"  results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.fasta | sed -e 's/_/ /g' | sort -nrk 6 | awk '$4>=2000 {print $0}' |sed -s 's/ /_/g' | sed -e 's/>//g' >  results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.5000-filtered-contigs.csv 
perl /home/swapnil/pipeline/tools/fastagrep.pl -f results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.5000-filtered-contigs.csv results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.fasta >  results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.fasta.filtered

awk '/^>/ {F = $1".fa"} {print > F}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.fasta.filtered
mv *.fa results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/
rename 's/>//' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/*.fa
#rename 's/_length.*/\.fa/' *.fa
ls results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/*.fa | sed "s/results\/25_blast_plasmid_finder_complete"$s"\/"$sp_p"\/raw_files\/$F1\/tmp\/fasta\///g" > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/list."$F1"_contigs.txt
##---------------------------
echo "sketching for $F1 contigs"
for contigs in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/list."$F1"_contigs.txt); do
(mash sketch -o results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/$contigs.msh results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/$contigs)> /dev/null 2>&1
done
##---------------------------
echo "Calculating mash distance for $F1"
for contigs in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/list."$F1"_contigs.txt); do
mash dist results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/$contigs.msh /home/swapnil/pipeline/tools/databases/plasmid_2/"$sp_p".msh | sed "s/\/1000//g" | awk '$3<0.2 && $5>200 {print "'$contigs'", $2, $3, $5}'| sort -k4,4rn | head -5 > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$contigs.1.matches.list

## if file is empty remove
if [ ! -s results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$contigs.1.matches.list ] ; then rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$contigs.1.matches.list
else
length=$(echo $contigs | awk -F'_' '{print $4}') #length
cov=$(echo $contigs | awk -F'_' '{print $6}' | sed 's/\.fa//g') #coverage
sed -i -e "s/$/ $length $cov/" results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$contigs.1.matches.list 
fi
done
##---------------------------
## take contigs with replicons (no criteria for mash distance or total number of mash)
if [ -s results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.WGS_replicon.contigs.list ] ; then
for contigs in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.WGS_replicon.contigs.list); do
mash dist results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/$contigs.msh /home/swapnil/pipeline/tools/databases/plasmid_2/"$sp_p".msh | sed "s/\/1000//g" | awk '{print "'$contigs'", $2, $3, $5}' | sort -k4,4rn | head -2 >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$contigs.2.matches.list
## if file is empty remove
if [ ! -s results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$contigs.2.matches.list ] ; then rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$contigs.2.matches.list
else
length=$(echo $contigs | awk -F'_' '{print $4}') #length
cov=$(echo $contigs | awk -F'_' '{print $6}' | sed 's/\.fa//g') #coverage
sed -i -e "s/$/ $length $cov/" results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$contigs.2.matches.list 
fi
done
fi

##---------------------------
Num_p_match=$(ls -U results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/ | grep -c matches.list)
if [ "$Num_p_match" -gt 0 ]; then
cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/*.matches.list | sort -u > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.csv

awk '{print $1}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.csv | sort -u > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.contig_names
awk '{print $2}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.csv | sort -u > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names

for plasmid in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names); do
awk '$2 == "'$plasmid'" {print $0}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.csv | awk '{sum+=$3; n++ } END {print "'$plasmid'", n, sum/n}' >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp1
awk '$2 == "'$plasmid'" {print $0}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.csv | awk -F'_' '{sum+=$4} END {print sum}' >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp2 #length
awk '$2 == "'$plasmid'" {print $0}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.csv | awk -F'_' '{sum+=$6; n++ } END {print sum/n}' >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp3 #coverage
done

paste results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp1 results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp2 results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp3 > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp4

cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp4 | sort -k2,2r > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.final.csv

awk '{print $1}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.final.csv > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names2

##---------------------------
mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs
for plasmid_contigs in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.contig_names); do
cp results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/$plasmid_contigs results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs/
done
cp -r results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmid_contigs
fi
##--------------------------
## Preliminary identification finished ----------------------------------------

if [ "$Num_p_match" -gt 0 ]; then
## stepwise plasmid identification round
(mkdir results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids) > /dev/null 2>&1
for multi_contigs in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names2) ; do
if [ "$(ls -A results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs/)" ] ; then 
awk '$2 == "'$multi_contigs'" {print $1}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.csv > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$multi_contigs.list
	for concat in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$multi_contigs.list); do
	if [ -f results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs/$concat ] ; then
	cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs/$concat >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$F1.$multi_contigs
	rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs/$concat
	fi
	done

	if [ -f results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$F1.$multi_contigs ] ; then
	A1=$(grep ">" results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$F1.$multi_contigs | tr "\n" " ")
	echo $multi_contigs $A1 >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$F1.contigs.list.csv
	sed -i 's/>.*/NNNNNNNNNN/g' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$F1.$multi_contigs
	sed -i "1i "'>'$F1.$multi_contigs"" results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$F1.$multi_contigs
	fi

	for concat in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$multi_contigs.list); do
	if [ -f results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs/$concat ] ; then
	(mash sketch -o results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/$F1.$multi_contigs.msh results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$F1.$multi_contigs)> /dev/null 2>&1
	mash dist results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta/$F1.$multi_contigs.msh /home/swapnil/pipeline/tools/databases/plasmid_2/*.msh | sed "s/\/1000//g" | awk '$3<0.2 && $5>200 {print "'$F1.$multi_contigs'", $2, $3, $5}'| sort -k4,4rn | head -5 > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$multi_contigs.matches.list2
	if [ ! -s results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$multi_contigs.matches.list2 ] ; then rm results/25_blast_plasmid_finder_completefall"$s"/$F1/tmp/$F1.$multi_contigs.matches.list2 
	else
	length=$(awk -F'_' '{sum+=$4} END {print sum}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$multi_contigs.list) #length
	cov=$(awk -F'_' '{sum+=$6} END {print sum}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$multi_contigs.list) #coverage
	sed -i -e "s/$/ $length $cov/" results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$multi_contigs.matches.list2
	fi
	fi	
	done
fi
done
fi
#===============================

## Stat

echo "Isolate Plasmid contigs_involved size(bp) Ref_size(bp) %covered pCov gCov replicon WGS_replicons contigs" > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.csv
## 
if [ "$Num_p_match" -gt 0 ]; then
ls results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/*.fa | awk -F'/' '{print $7}' | sed 's/\.fa//g' > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/final_plasmid_names.list

for final_plasmid in $(cat results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/final_plasmid_names.list); do
if [ ! -z "$final_plasmid" ]; then
name=$(echo $final_plasmid | sed 's/'$F1'\.//g' | sed 's/\.fa//g' )
## find replicon
if [ -f results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$final_plasmid.fa ] ; then
blastn -db /home/swapnil/pipeline/tools/databases/plasmid_replicons/310519_all.fsa -query results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$final_plasmid.fa -out results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/plasmid_replicon/$final_plasmid.blast.tmp -max_target_seqs 1 -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
replicon=$(grep $final_plasmid results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/plasmid_replicon/$final_plasmid.blast.tmp | awk '{print $1}' | sed '/^$/d' | tr "\n" "_" )
## Number of contigs forming reconstructed_plasmid
Num_cont=$(grep -c "NNNNNNNNNN" results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$final_plasmid.fa)
## size of the reconstructed_plasmid
p_size=$(bioawk -c fastx '{ print $name, length($seq) }' <results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$final_plasmid.fa | awk 'NR==1 {print $2}' )
## contigs constructing plasmid
CCP=$(grep $name results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/final_plasmids/$F1.contigs.list.csv | awk '{ s = ""; for (i = 2; i <= NF; i++) s = s $i " "; print s }' | sed "s/>/ /g")
## coverage of the reconstructed_plasmid 
cov=$(awk -F'_' '{sum+=$6} END {print sum}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.$name.fa.list)
## size of the closest plasmid
V2=$(grep $name /home/swapnil/pipeline/tools/databases/plasmid_2/00_All_plasmids_v16012019.fna.list | awk 'NR==1 {print $2}') 
## 100*(size of reconstructed_plasmid / size of closest plasmid)
V3=$(echo "scale=2; 100*$p_size/$V2" | bc)
fi
fi

## WGS replicon
WGS_replicon=$(grep $F1 results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.WGS_replicon.list | awk '{print $2}')

if [ -z "$name" ]; then name=$(echo "NA") ; fi
if [ -z "$replicon" ]; then replicon=$(echo "NA") ; fi
if [ -z "$Num_cont" ]; then Num_cont=$(echo "NA") ; fi
if [ -z "$p_size" ]; then p_size=$(echo "NA") ; fi
if [ -z "$CCP" ]; then CCP=$(echo "NA") ; fi
if [ -z "$cov" ]; then cov=$(echo "NA") ; fi
if [ -z "$V2" ]; then V2=$(echo "NA") ; fi
if [ -z "$V3" ]; then V3=$(echo "NA") ; fi
if [ -z "$WGS_replicon" ]; then WGS_replicon=$(echo "NA") ; fi

echo $F1 $name $Num_cont $p_size $V2 $V3 $cov $genome_cov $replicon >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.1.csv
echo $CCP >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.2.csv
done
else
echo "$F1 NA NA NA NA NA NA NA NA" > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.1.csv
echo "NA" > results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.2.csv
WGS_replicon=$(grep $F1 results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/$F1.WGS_replicon.list | awk '{print $2}')
echo $WGS_replicon >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.3.csv
fi



paste results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.1.csv results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.3.csv results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.2.csv >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.csv

awk 'FNR>1 {print $0}' results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/$F1.results.csv >> results/25_blast_plasmid_finder_complete"$s"/$sp_p/results.csv

#===============================
## remove unncessary files
#rm -rf results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/fasta
#rm -rf results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/plasmid_contigs
#rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/*.list
#rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp1
#rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp2
#rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp3
#rm results/25_blast_plasmid_finder_complete"$s"/$sp_p/raw_files/$F1/tmp/all.matches.list.plasmid_names.averaged.tmp4
#===============================

done

ssconvert results/25_blast_plasmid_finder_complete"$s"/$sp_p/results.csv results/25_blast_plasmid_finder_complete"$s"/$sp_p/results.csv.xlsx
###############################################################################
echo "Completed.. plasmid serach complete ------------------------------------"

###############################################################################
