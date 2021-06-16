#!/bin/bash
###############################################################################
#7 map reads to reference
###############################################################################
echo "started.... step-7 mapping reads to the reference -----------------------"
###############################################################################
## preliminary file preparations

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | sed 's/list\.//g' | sed 's/\.txt//g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany)"
    read path_tmp
    if [ -z $path_tmp ] ; then
        path=$(echo "$path_tmp")
        else
        path=$(echo "$path_tmp/")
    fi

    (mkdir results/10_read_mapping) > /dev/null 2>&1
    (mkdir results/10_read_mapping/tmp) > /dev/null 2>&1

for F1 in $(cat $list) ; do
    echo "running.... step-7 mapping reads to the reference for..." $F1
    mkdir results/10_read_mapping/$F1
    bwa index results/00_ref/ref.*.fasta  
    bwa mem results/00_ref/ref.*.fasta "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz  "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz > results/10_read_mapping/$F1/$F1.sam 
    samtools view -bS  results/10_read_mapping/$F1/$F1.sam  > results/10_read_mapping/$F1/$F1.bam  
    samtools sort results/10_read_mapping/$F1/$F1.bam > results/10_read_mapping/$F1/"$F1"_sorted.bam
    samtools flagstat results/10_read_mapping/$F1/"$F1"_sorted.bam > results/10_read_mapping/tmp/$F1.7_map_reads_to_reference.statistics.tab
    echo "finished... step-7 mapping reads to the reference for..." $F1
done

###############################################################################
echo "completed.... step-7 mapping reads to the reference ---------------------"
###############################################################################

samtools mpileup -uf results/00_ref/ref.*.fasta results/10_read_mapping/$F1/"$F1"_sorted.bam | bcftools view -bvcg - > results/10_read_mapping/$F1/var.$F1.raw.bcf  
bcftools view results/10_read_mapping/$F1/var.$F1.raw.bcf | vcfutils.pl varFilter -D30 > results/10_read_mapping/$F1/var.$F1.flt.vcf
gzip -c results/10_read_mapping/$F1/var.$F1.flt.vcf > results/10_read_mapping/$F1/var.$F1.flt.vcf.gz
bcftools consensus -f results/00_ref/ref.*.fasta results/10_read_mapping/$F1/var.$F1.flt.vcf.gz > results/10_read_mapping/$F1/$F1.consensus.fa
# or
samtools mpileup -uf results/00_ref/ref.*.fasta results/10_read_mapping/SE/SE_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > results/10_read_mapping/SE/SE.cns.fq
seqret -osformat fasta results/10_read_mapping/$F1/$F1.cns.fq -out2 results/10_read_mapping/$F1/$F1.cns.fq.consensus.fa

#or 
samtools mpileup -uf results/00_ref/ref.*.fasta results/10_read_mapping/$F1/"$F1"_sorted.bam | bcftools call -mv -Oz -o results/10_read_mapping/$F1/$F1.calls.vcf.gz