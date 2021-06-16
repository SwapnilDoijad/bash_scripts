#!/bin/bash
###############################################################################

echo "started to arrange closed genome as per dnaA gene -----------------------"

###############################################################################
    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate myenv

echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

#------------------------------------------------------------------------------v
(mkdir results/04_assembly)> /dev/null 2>&1
(mkdir results/04_assembly/dnaA_arranged)> /dev/null 2>&1
#------------------------------------------------------------------------------
for F1 in $(cat $list); do
    (mkdir results/04_assembly/dnaA_arranged/$F1)> /dev/null 2>&1
    (mkdir results/04_assembly/dnaA_arranged/$F1/tmp)> /dev/null 2>&1
    (mkdir results/04_assembly/dnaA_arranged/$F1/fasta_seperated)> /dev/null 2>&1
    sed -i 's/>.*/>'$F1'/g' results/04_assembly/raw_files/$F1/assembly.fasta 
    cp results/04_assembly/raw_files/$F1/assembly.fasta results/04_assembly/dnaA_arranged/$F1/tmp/$F1.unarranged.fasta
    sed -i "s/ .*//g" results/04_assembly/dnaA_arranged/$F1/tmp/$F1.unarranged.fasta
    bioawk -c fastx '{ print $name, length($seq) }' results/04_assembly/dnaA_arranged/$F1/tmp/$F1.unarranged.fasta | sort -nrk 2 | awk '{print $1}' > results/04_assembly/dnaA_arranged/$F1/tmp/$F1.contigs.list
    chr=$(bioawk -c fastx '{ print $name, length($seq) }' results/04_assembly/dnaA_arranged/$F1/tmp/$F1.unarranged.fasta | sort -nrk 2 | head -1 | awk '{print $1}')
    perl /home/swapnil/pipeline/tools/split_multifasta.pl --in results/04_assembly/dnaA_arranged/$F1/tmp/$F1.unarranged.fasta --output_dir=results/04_assembly/dnaA_arranged/$F1/fasta_seperated   
    fastahack -i results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa
    echo "running annotation for $F1 $chr"
    prokka --quiet --force --outdir results/04_assembly/dnaA_arranged/$F1/tmp --prefix $F1 results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa
    echo "finished annotation for $F1 $chr"
    checkStrand=$(grep "dnaA" results/04_assembly/dnaA_arranged/$F1/tmp/$F1.gff | awk 'NR==1 {print $7}')

    if [ "$checkStrand" == "+" ]; then
        V1=$(grep "dnaA" results/04_assembly/dnaA_arranged/$F1/tmp/$F1.gff | awk 'NR==1 {print $4}')
        if [ "$V1" -gt 100 ] ; then
            V2=$(( $V1 - 100 ))
            V3=$(( $V1 - 101 ))
            V4=$(bioawk -c fastx '{ print length($seq) }' <results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa)
            fastahack -r "$chr":"$V2".."$V4" results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa > results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-end.tmp
            fastahack -r "$chr":1.."$V3" results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa > results/04_assembly/dnaA_arranged/$F1/tmp/$F1.start-dnaA.tmp
            cat results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-end.tmp results/04_assembly/dnaA_arranged/$F1/tmp/$F1.start-dnaA.tmp > results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-arranged.fasta
            sed -i "1i "'>'$chr"" results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-arranged.fasta
            cp results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-arranged.fasta results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa
            else
            cp results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-arranged.fasta
            cp results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-arranged.fasta results/04_assembly/dnaA_arranged/$F1.chr.fasta
        fi
        else
        echo "dnaA is on negative strand for $F1"
        X1=$(grep "dnaA" results/04_assembly/dnaA_arranged/$F1/tmp/$F1.gff | awk 'NR==1 {print $5}')
        X2=$(( $X1 + 100 ))
        X3=$(( $X1 + 101 ))
        X4=$(bioawk -c fastx '{ print length($seq) }' <results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa)
        fastahack -r "$chr":1.."$X2" results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa > results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-end.-RC.tmp
        fastahack -r "$chr":"$X3".."$X4" results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa > results/04_assembly/dnaA_arranged/$F1/tmp/$F1.end-dnaA.-RC.tmp
        cat results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-end.-RC.tmp results/04_assembly/dnaA_arranged/$F1/tmp/$F1.end-dnaA.-RC.tmp | perl -pe 'chomp;tr/ACGTNacgtn/TGCANtgcan/;$_=reverse."\n"' > results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-arranged.fasta
        sed -i "1i "'>'$chr"" results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-arranged.fasta
        cp results/04_assembly/dnaA_arranged/$F1/tmp/$F1.dnaA-arranged.fasta results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$chr.fsa
    fi

    for F2 in $(cat results/04_assembly/dnaA_arranged/$F1/tmp/$F1.contigs.list); do
    cat results/04_assembly/dnaA_arranged/$F1/fasta_seperated/$F2.fsa >> results/04_assembly/dnaA_arranged/$F1/$F1.fasta
    done
    sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/dnaA_arranged/$F1/$F1.fasta
    sed -i "1i "'>'$F1"" results/04_assembly/dnaA_arranged/$F1/$F1.fasta 
    sed -i '2d' results/04_assembly/dnaA_arranged/$F1/$F1.fasta
done

conda deactivate
###############################################################################

echo "started to arrange closed genome as per dnaA gene -----------------------"

###############################################################################
