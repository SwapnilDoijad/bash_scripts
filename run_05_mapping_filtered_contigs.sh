#!/bin/bash
###############################################################################
#6 map filtered contigs
###############################################################################
echo "started...... step-6 mapping filtered contigs to the reference ---------"
###############################################################################
    alias python=python3
    (rm /home/swapnil/tools/mauve/*.fasta) > /dev/null 2>&1

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | sed 's/list\.//g' | sed 's/\.txt//g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "method? for medusa type m for mauve type ma"
    read method

    (mkdir tmp) > /dev/null 2>&1
    (mkdir results/00_ref) > /dev/null 2>&1
    (mkdir results/00_ref/00_mash_sketch) > /dev/null 2>&1
    (mkdir results/04_assembly/all_fasta/mash_sketch) > /dev/null 2>&1
    (mkdir data/references/medusa) > /dev/null 2>&1

     WorDir=$(echo $PWD)
###############################################################################
## choose closest reference
echo "preparing closest reference"
ls data/references/*.gbk | awk -F '/' '{print $NF}' | sed 's/\.gbk//g' > list.ref.txt
for ref in $(cat list.ref.txt); do
/home/swapnil/pipeline/tools/GB2Fasta.pl data/references/$ref.gbk results/00_ref/$ref.fasta
sed -i '1i >'$ref'' results/00_ref/$ref.fasta
(mash sketch -o results/00_ref/00_mash_sketch/$ref.msh results/00_ref/$ref.fasta) > /dev/null 2>&1 
done

for F1 in $(cat $list); do
(mash sketch -o results/04_assembly/all_fasta/mash_sketch/$F1.msh results/04_assembly/all_fasta/$F1.fasta) > /dev/null 2>&1 
done

for F1 in $(cat $list);do
(mash dist results/04_assembly/all_fasta/mash_sketch/$F1.msh results/00_ref/00_mash_sketch/*.msh | sort -k 3,3 > results/04_assembly/all_fasta/mash_sketch/$F1.msh.dist.csv) > /dev/null 2>&1 
done
###############################################################################
for F1 in $(cat $list); do
    if [ "$method" == "ma" ] ; then

        echo "running.... step-6 mapping filtered contigs to the reference by mauve for..." $F1

        ##choose reference
        new_ref=$(head -1 results/04_assembly/all_fasta/mash_sketch/$F1.msh.dist.csv | awk '{print $2}' | awk -F '/' '{print $NF}' | sed 's/\.fasta//g')
        echo "reference $new_ref"
        cp data/references/$new_ref.gbk data/references/medusa/$new_ref.gbk

        cp results/04_assembly/raw_files/$F1/$F1.fasta /home/swapnil/tools/mauve/

        WorDir=$(echo $PWD)
        cd /home/swapnil/tools/mauve
        (java -Xmx5000m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output $F1 -ref $WorDir/results/00_ref/$new_ref.fasta -draft $F1.fasta) > /dev/null 2>&1 
        cd $WorDir

        Var1=$(ls /home/swapnil/tools/mauve/$F1/ | sort -r | head -1)
        cp /home/swapnil/tools/mauve/$F1/$Var1/$F1.fasta  results/04_assembly/raw_files/$F1/"$F1".aligned.fasta
        cp results/04_assembly/raw_files/$F1/"$F1".aligned.fasta results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta
        sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta
        sed -i "1i "'>'$F1"" results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta

        (rm -rf /home/swapnil/tools/mauve/$F1) > /dev/null 2>&1
        (rm /home/swapnil/tools/mauve/$F1.fasta) > /dev/null 2>&1

    elif [ "$method" == "m" ] ; then
        echo "running.... step-6 mapping filtered contigs to the reference by medusa for..." $F1

        if [ -f results/04_assembly/raw_files/$F1/$F1.fasta ] ; then
            cp results/04_assembly/raw_files/$F1/$F1.fasta /home/swapnil/tools/medusa/
            else        
            cp results/04_assembly/raw_files/$F1/contigs.fasta /home/swapnil/tools/medusa/$F1.contigs.fasta
        fi

        new_ref=$(head -1 results/04_assembly/all_fasta/mash_sketch/$F1.msh.dist.csv | awk '{print $2}' | awk -F '/' '{print $NF}' | sed 's/\.fasta//g')
        echo "reference $new_ref"
        cp data/references/$new_ref.gbk data/references/medusa/$new_ref.gbk

        cd /home/swapnil/tools/medusa
        (java -jar medusa.jar -f $WorDir/data/references/medusa -i $F1.*.fasta -d 10 ) > /dev/null 2>&1 
        cd $WorDir

        ##sometime  medus fail, so if medusa fail then run with mauve
        if [ -f /home/swapnil/tools/medusa/"$F1".*.fastaScaffold.fasta ] ; then
            mv /home/swapnil/tools/medusa/"$F1".*.fastaScaffold.fasta results/04_assembly/raw_files/$F1/
            ( rm /home/swapnil/tools/medusa/*.delta ) > /dev/null 2>&1
            ( rm /home/swapnil/tools/medusa/*.coords ) > /dev/null 2>&1
            ( rm /home/swapnil/tools/medusa/*.ntref ) > /dev/null 2>&1
            ( rm /home/swapnil/tools/medusa/*.mgaps ) > /dev/null 2>&1
            cp results/04_assembly/raw_files/$F1/"$F1".*.fastaScaffold.fasta results/04_assembly/raw_files/$F1/"$F1".aligned.fasta
            cp results/04_assembly/raw_files/$F1/"$F1".aligned.fasta results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta
            sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta
            sed -i "1i "'>'$F1"" results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta

            mv /home/swapnil/tools/medusa/$F1.*.fasta_SUMMARY results/04_assembly/raw_files/$F1/
            rm /home/swapnil/tools/medusa/$F1.*
            rm data/references/medusa/$new_ref.gbk
            else
            echo "Medusa failed, running.... step-6 mapping filtered contigs to the reference by mauve for... $F1"

                ##choose reference
                new_ref=$(head -1 results/04_assembly/all_fasta/mash_sketch/$F1.msh.dist.csv | awk '{print $2}' | awk -F '/' '{print $NF}' | sed 's/\.fasta//g')
                echo "reference $new_ref"
                cp data/references/$new_ref.gbk data/references/medusa/$new_ref.gbk

                cp results/04_assembly/raw_files/$F1/$F1.fasta /home/swapnil/tools/mauve/

                WorDir=$(echo $PWD)
                cd /home/swapnil/tools/mauve
                (java -Xmx5000m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output $F1 -ref $WorDir/results/00_ref/$new_ref.fasta -draft $F1.fasta) > /dev/null 2>&1 
                cd $WorDir

                Var1=$(ls /home/swapnil/tools/mauve/$F1/ | sort -r | head -1)
                cp /home/swapnil/tools/mauve/$F1/$Var1/$F1.fasta  results/04_assembly/raw_files/$F1/"$F1".aligned.fasta
                cp results/04_assembly/raw_files/$F1/"$F1".aligned.fasta results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta
                sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta
                sed -i "1i "'>'$F1"" results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta

                (rm -rf /home/swapnil/tools/mauve/$F1) > /dev/null 2>&1
                (rm /home/swapnil/tools/mauve/*.fasta) > /dev/null 2>&1

        fi
    fi

        (mkdir results/04_assembly/all_fasta) > /dev/null 2>&1


        ## reformat fasta 
        fasta_formatter -i results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta -w 60 -o results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta.tmp
        mv results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta.tmp results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta

        cp results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta results/04_assembly/all_fasta/$F1.fasta

        if [ ! -f results/04_assembly/all_fasta/$F1.fasta ]; then
        echo "failed! mapping contigs for $F1.fasta" 
        echo "failed! mapping contigs for $F1.fasta" >> results/failed_list.txt
        fi

done
###############################################################################
echo "completed.... step-6 mapping filtered contigs to the reference ---------"
###############################################################################
