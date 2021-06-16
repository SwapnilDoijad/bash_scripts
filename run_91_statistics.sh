#!/bin/bash
###############################################################################
#10 statistics

echo "started step-10 statistics ----------------------------------------------"

###############################################################################
rm results/statistics.csv

if [ ! -f results/statistics.csv ]; then
touch results/statistics.csv
fi


if grep -Fxq  "strain raw-F-reads raw-R-reads total-raw-read average-read-length-raw-reads filtered-F-reads filtered-R-reads Average-read-length-filtered-reads filtered-total-reads(%-of-total-raw-reads) mapped-reads(%-of-total-filtered-reads)(%-of-total-raw-reads) total-(unfiltered)-contigs final-filtered-(>500bp>10x)-contigs CoverageBasedOnAssembledContigs GenomeSize(bp) CoverageBasedOnAssembledGenome CoverageBasedOnRefSize %SizeToRefGenomeSize GC% rRNA tRNA genes CDS SNPs-to-Ref(Inter+Intra) dN/dS ts/tv" results/statistics.csv; then
:
else
ex -sc '1i|strain raw-F-reads raw-R-reads total-raw-read average-read-length-raw-reads filtered-F-reads filtered-R-reads Average-read-length-filtered-reads filtered-total-reads(%-of-total-raw-reads) mapped-reads(%-of-total-filtered-reads)(%-of-total-raw-reads) total-(unfiltered)-contigs final-filtered-(>500bp>10x)-contigs CoverageBasedOnAssembledContigs GenomeSize(bp) CoverageBasedOnAssembledGenome CoverageBasedOnRefSize %SizeToRefGenomeSize GC% rRNA tRNA genes CDS SNPs-to-Ref(Inter+Intra) dN/dS ts/tv' -cx results/statistics.csv
fi

###############################################################################

for F1 in $(cat list.txt)
do
V1=$(awk 'FNR == 1 {print $1}' results/01_raw_read_count/raw_files/$F1/tmp/$F1.1_raw_read_count.statistics.tab)			#raw-F-reads
V2=$(awk 'FNR == 1 {print $2}' results/01_raw_read_count/raw_files/$F1/tmp/$F1.1_raw_read_count.statistics.tab)			#raw-R-reads
V3=$(expr $V1 + $V2 )														#total-read
V3a=$(awk 'FNR == 1 {print $3}' results/01_raw_read_count/raw_files/$F1/tmp/$F1.1_raw_read_count.statistics.tab)			#Average-read-length-raw-reads

V4=$(awk 'FNR == 1 {print $1}' results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.statistics.tab)		#filtered-F-reads
V5=$(awk 'FNR == 1 {print $2}' results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.statistics.tab)		#filtered-R-reads
V6a=$(expr $V4 + $V5 )									#filtered-total-reads
V5a=$(awk 'FNR == 1 {print $3}' results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.statistics.tab)		#Average-read-length-filtered-reads								
V6b=$((100 * $V6a / $V3))
V6c=$V6a'('$V6b%')'										#filtered-total-reads & %-of-total-raw-reads
#V7a=$(awk 'FNR == 5 {print $1}' results/tmp/$F1.7_map_reads_to_reference.statistics.tab)	
#V7b=$(( 100 * $V7a / $V6a ))
#V7c=$(( 100 * $V7a / $V3 ))
#V7d=$V7a'('$V7b%')''('$V7c%')'                 						#mapped-reads(%-of-total-filtered-reads)(%-of-total-raw-reads)
V8=$(awk 'FNR == 1 {print $1}' results/04_assembly/raw_files/$F1/tmp/$F1.5_filtering_contigs.statistics.tab)		#total-(unfiltered)-contig
V9=$(awk 'FNR == 1 {print $2}' results/04_assembly/raw_files/$F1/tmp/$F1.5_filtering_contigs.statistics.tab)		#final-filtered-(>500bp>10x)-contigs
V10=$(awk 'FNR == 1 {print $3}' results/04_assembly/raw_files/$F1/tmp/$F1.5_filtering_contigs.statistics.tab)		#ContigCoverage
	
###############################################################################

G1=$(awk 'FNR == 1 {print $2}' results/08_annotation/raw_files/$F1/tmp/$F1.8_annotation.statistics.txt) 			#GenomeSize(bp)
G2=$(awk 'FNR == 2 {print $2}' results/08_annotation/raw_files/$F1/tmp/$F1.8_annotation.statistics.txt)			#rRNA
G3=$(awk 'FNR == 3 {print $2}' results/08_annotation/raw_files/$F1/tmp/$F1.8_annotation.statistics.txt)			#tRNA
G4=$(awk 'FNR == 4 {print $2}' results/08_annotation/raw_files/$F1/tmp/$F1.8_annotation.statistics.txt)			#genes
G5=$(awk 'FNR == 5 {print $2}' results/08_annotation/raw_files/$F1/tmp/$F1.8_annotation.statistics.txt)			#cds
G6=$(awk "BEGIN {print ($V6a*$V5a/$G1); exit}")							#GenomeCoverage = filtered-total-reads * Average-read-length-filtered-reads / GenomeSize(bp)
G6b=$(awk 'FNR == 6 {print $2}' results/08_annotation/raw_files/$F1/tmp/$F1.8_annotation.statistics.txt)
###############################################################################

#G7a=$(wc -c < results/00_ref/ref.*.fasta)							#
#G7b=$(awk "BEGIN {print ($V6a*$V5a/$G7a); exit}")
#G7c=$(awk "BEGIN {print ($G1*100/$G7a); exit}" | awk '{printf "%.2f\n", $1, $2}')		#%SizeToRefGenomeSize		
#G8=$(awk 'FNR == 1 {print $2}' results/tmp/$F1.tstv.txt) 					#ts/tv
###############################################################################

#G9=$(awk 'FNR == 1 {print $1}' results/tmp/$F1.13_snippy.SNPs-statistics.txt)			#SNPs-to-Ref(Inter+Intra)
#G10=$(awk 'FNR == 1 {print $2}' results/tmp/$F1.13_snippy.SNPs-statistics.txt)

###############################################################################

echo $F1 $V1 $V2 $V3 $V3a $V4 $V5 $V5a $V6c $V8 $V9 $V10 $G1 $G6 $G7b >> results/statistics.csv #$G7c $G6b $G2 $G3 $G4 $G5 $G8 

###############################################################################

done

sed -i 's/ /\t/g' results/statistics.csv

unoconv -i FilterOptions=09,,system,1 -f xls results/statistics.csv

echo "completed step-10 statistics --------------------------------------------"

###############################################################################








