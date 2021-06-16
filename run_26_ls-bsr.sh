#!/bin/bash
###############################################################################
#26 ls-bsr

###############################################################################

echo "started.... ls-bsr -----------------------------------------------------"

###############################################################################
## preliminary file preparation

echo "provide list file (for e.g. all)"
echo "---------------------------------------------------------------------"
ls list.*.txt | sed 's/ /\n/g'
echo "---------------------------------------------------------------------"
read l
list=$(echo "list.$l.txt")

echo "provide complete path to genes.fasta (for e.g. data/genes_to_blast/geneA.fasta)"
read geneFasta_path
sfx_tmp=$(echo $geneFasta_path | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}')
sfx=$(echo "_$sfx_tmp")

if [ ! -d results/26_ls-bsr$sfx ]; then
mkdir results/26_ls-bsr$sfx
fi

if [ ! -d results/26_ls-bsr$sfx/fasta ]; then
mkdir results/26_ls-bsr$sfx/fasta
fi

for F1 in $(cat $list); do
cp results/04_assembly/all_fasta/"$F1".fasta results/26_ls-bsr$sfx/fasta/
done

###############################################################################

prefix=$(echo $geneFasta_path | awk -F "/" '{print $NF}' | awk -F"." '{print $1}')

python /home/swapnil/tools/ls_bsr/ls_bsr.py -d results/26_ls-bsr$sfx/fasta/ -g $geneFasta_path -y results/26_ls-bsr$sfx/tmp -x $prefix -b blastn

mv $prefix* results/26_ls-bsr$sfx

awk '{for (i=1;i<=NF;i++) sum[i]+=$i;}; END{for (i in sum) print "for column "i" is " sum[i];}' results/26_ls-bsr$sfx/"$prefix"_bsr_matrix.txt | sort -k 5,5 > results/26_ls-bsr$sfx/"$prefix"_bsr_matrix.summed-sorted.txt 

awk '{print $5}' results/26_ls-bsr$sfx/"$prefix"_bsr_matrix.summed-sorted.txt > results/26_ls-bsr$sfx/tmp1.tmp
cat results/26_ls-bsr$sfx/tmp1.tmp | uniq > results/26_ls-bsr$sfx/tmp2.tmp

for F1 in $(cat results/26_ls-bsr$sfx/tmp2.tmp); do
V1=$(grep -cw "$F1" results/26_ls-bsr$sfx/tmp1.tmp);
echo $F1 $V1 >> results/26_ls-bsr$sfx/"$prefix"_bsr_matrix.summed-sorted.summary.txt 
done

rm -rf results/26_ls-bsr$sfx/fasta/

###############################################################################

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' results/26_ls-bsr$sfx/"$prefix"_bsr_matrix.txt > results/26_ls-bsr$sfx/"$prefix"_bsr_matrix.transposed.txt

###############################################################################

echo "Completed.... ls-bsr ---------------------------------------------------"

###############################################################################
