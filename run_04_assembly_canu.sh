#!/bin/bash
###############################################################################
#4-B long read assembly

###############################################################################

echo "started.... step-4b long read assembly ---------------------------------"

###############################################################################

for F1 in $(cat list.txt); do

if [ ! -d results/04_assembly_long-reads_canu ]; then
mkdir results/04_assembly_long-reads_canu
fi

perl /home/swapnil/tools/canu/Linux-amd64/bin/canu -p $F1 -d results/04_assembly_long-reads_canu/$F1 genomesize=5.0m -pacbio-raw data/pacbio/$F1/"$F1"_1.fastq

if [ ! -d results/04_assembly_long-reads_canu/tmp ]; then
mkdir results/04_assembly_long-reads_canu/tmp
fi
#------------------------------------------------------------------------------
## arranging the contigs from higher length to lower length
cp results/04_assembly_long-reads_canu/$F1/$F1.contigs.fasta results/04_assembly_long-reads_canu/$F1/$F1.contigs.2.fasta
sed -i 's/ /_/g' results/04_assembly_long-reads_canu/$F1/$F1.contigs.2.fasta
grep ">" results/04_assembly_long-reads_canu/$F1/$F1.contigs.2.fasta | sort -nrk 2 -t'=' > results/04_assembly_long-reads_canu/tmp/$F1.contigs-co-ordinates.arranged.txt
for F2 in $(cat results/04_assembly_long-reads_canu/tmp/$F1.contigs-co-ordinates.arranged.txt); do
V1=$(echo $F2 | sed 's/>//g')
perl /home/swapnil/pipeline/tools/fastagrep.pl $V1 results/04_assembly_long-reads_canu/$F1/$F1.contigs.2.fasta >> results/04_assembly_long-reads_canu/$F1/$F1.contigs.arranged.fasta 
done
#------------------------------------------------------------------------------
## mapping the reads using bwa for pilon
mkdir results/04_assembly_long-reads_canu/$F1/read_mapping
bwa index results/04_assembly_long-reads_canu/$F1/$F1.contigs.arranged.fasta 
bwa mem results/04_assembly_long-reads_canu/$F1/$F1.contigs.arranged.fasta results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz > results/04_assembly_long-reads_canu/$F1/read_mapping/$F1.sam 
samtools view -bS  results/04_assembly_long-reads_canu/$F1/read_mapping/$F1.sam  > results/04_assembly_long-reads_canu/$F1/read_mapping/$F1.bam  
samtools sort results/04_assembly_long-reads_canu/$F1/read_mapping/$F1.bam > results/04_assembly_long-reads_canu/$F1/read_mapping/"$F1"_sorted.bam
samtools flagstat results/04_assembly_long-reads_canu/$F1/read_mapping/"$F1"_sorted.bam > results/04_assembly_long-reads_canu/$F1/read_mapping/$F1.7_map_reads_to_reference.statistics.tab
#------------------------------------------------------------------------------
## running pilon-for-polishing 
samtools index results/04_assembly_long-reads_canu/$F1/read_mapping/"$F1"_sorted.bam
java -jar /home/swapnil/tools/pilon-1.22.jar --genome results/04_assembly_long-reads_canu/$F1/$F1.contigs.arranged.fasta --frags results/04_assembly_long-reads_canu/$F1/read_mapping/"$F1"_sorted.bam --output $F1 --outdir results/04_assembly_long-reads_canu/$F1/pilon --changes --vcf
#------------------------------------------------------------------------------
cp results/04_assembly_long-reads_canu/$F1/pilon/$F1.fasta results/04_assembly_long-reads_canu/$F1/$F1.joined.fasta
sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly_long-reads_canu/$F1/$F1.joined.fasta
sed -i "1i "'>'$F1"" results/04_assembly_long-reads_canu/$F1/$F1.joined.fasta
sed -i "\$aNNNNNNNNNN" results/04_assembly_long-reads_canu/$F1/$F1.joined.fasta 
 
done
###############################################################################

echo "completed.... step-4b long read assembly -------------------------------"

###############################################################################
