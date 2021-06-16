#!/bin/bash
###############################################################################
echo "started.... ANI-enveomics ----------------------------------------------"
###############################################################################
## initial inout and file preparation
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")
    sfx=$(echo "_$l")

    echo "Type N for ANI or type A for AAI" 
    read answer

    (mkdir results/37_ANI-enveomics"$sfx") > /dev/null 2>&1
    (mkdir results/37_ANI-enveomics"$sfx"/tmp) > /dev/null 2>&1

###############################################################################
if [ "$answer" == "N" ]; then
    ## running ANI
    for F1 in $(cat $list); do
            (mkdir results/37_ANI-enveomics"$sfx"/$F1) > /dev/null 2>&1
        echo "subj quer two-way-ANI(%) TWA-fragm subj-fragm quer-fragm %Subj_covered %Quer_covered" > results/37_ANI-enveomics"$sfx"/$F1.ANI.csv

        for F2 in $(cat $list); do
        echo "running $F1 vs $F2"
            ( ruby /home/swapnil/tools/enveomics-master/Scripts/ani.rb --seq1 results/04_assembly/all_fasta/$F1.fasta --seq2 results/04_assembly/all_fasta/$F2.fasta --res results/37_ANI-enveomics"$sfx"/$F1/$F1.$F2.res.out_file --quiet --threads 8 ) > /dev/null 2>&1

            V1=$(awk 'FNR==3{print $3}' results/37_ANI-enveomics"$sfx"/$F1/$F1.$F2.res.out_file | sed 's/\%//g')
            V2=$(awk 'FNR==3{print $7}' results/37_ANI-enveomics"$sfx"/$F1/$F1.$F2.res.out_file)
            V3=$(awk 'FNR==1{print $8}' results/37_ANI-enveomics"$sfx"/$F1/$F1.$F2.res.out_file)
            V4=$(awk 'FNR==2{print $8}' results/37_ANI-enveomics"$sfx"/$F1/$F1.$F2.res.out_file)
            V5=$(echo 'scale=1;100*'$V2'/'$V3'' | bc )
            V6=$(echo 'scale=1;100*'$V2'/'$V4'' | bc )           
            echo "$F1 $F2 $V1 $V2 $V3 $V4 $V5 $V6" >> results/37_ANI-enveomics"$sfx"/$F1.ANI.csv
        done
    done

    ##-------------------------------------------------------------------------
    echo "generating matrix"
    ## get first columum, use of $F1 is OK
    cat results/37_ANI-enveomics"$sfx"/$F1.ANI.csv | awk '{print $2}' > results/37_ANI-enveomics"$sfx"/tmp/00.second_column_isolate_list.tmp
    cat results/37_ANI-enveomics"$sfx"/$F1.ANI.csv | awk 'NR>1 {print $2}' > results/37_ANI-enveomics"$sfx"/tmp/00.second_column_isolate_list.tmp.list
    cp results/37_ANI-enveomics"$sfx"/tmp/00.second_column_isolate_list.tmp results/37_ANI-enveomics"$sfx"/tmp/matrix.tmp
    for F1 in $(cat results/37_ANI-enveomics"$sfx"/tmp/00.second_column_isolate_list.tmp.list); do
    cat results/37_ANI-enveomics"$sfx"/$F1.ANI.csv | awk 'NR>1 {print $3}' | sed '1 i\'$F1'' > results/37_ANI-enveomics"$sfx"/tmp/$F1.tmp
    paste results/37_ANI-enveomics"$sfx"/tmp/matrix.tmp results/37_ANI-enveomics"$sfx"/tmp/$F1.tmp > results/37_ANI-enveomics"$sfx"/tmp/matrix.tmp2
    cp results/37_ANI-enveomics"$sfx"/tmp/matrix.tmp2 results/37_ANI-enveomics"$sfx"/tmp/matrix.tmp
    done
    sed -e 's/quer//g' results/37_ANI-enveomics"$sfx"/tmp/matrix.tmp > results/37_ANI-enveomics"$sfx"/matrix.tab
    ##-------------------------------------------------------------------------
    ## For R
    
    sed -i "s/37_ANI-enveomics/37_ANI-enveomics"$sfx"/g" /home/swapnil/pipeline/tools/plot_and_tree.TWA-distance.r
    Rscript /home/swapnil/pipeline/tools/plot_and_tree.TWA-distance.r
    sed -i "s/37_ANI-enveomics"$sfx"/37_ANI-enveomics/g" /home/swapnil/pipeline/tools/plot_and_tree.TWA-distance.r
    (mv Rplots.pdf results/37_ANI-enveomics"$sfx"/)> /dev/null 2>&1
    (mv tree.nwk results/37_ANI-enveomics"$sfx"/)> /dev/null 2>&1



    else
    ## runnign AAI
    for F1 in $(cat $list); do
            (mkdir results/37_AAI-enveomics"$sfx") > /dev/null 2>&1
        echo "subj quer two-way-AAI TWA-fragm subj-fragm quer-fragm" > results/37_AAI-enveomics"$sfx".AAI.csv

        for F2 in $(cat $list); do
        echo "running $F1 vs $F2"
            (ruby /home/swapnil/tools/enveomics-master/Scripts/aai.rb --seq1 results/08_annotation/all_faa/$F1.faa --seq2 results/08_annotation/all_faa/$F2.faa --res results/37_AAI-enveomics"$sfx"/$F1.$F2.res.out_file --rbm results/37_AAI-enveomics"$sfx"/$F1.$F2.rbm.tab --quiet) > /dev/null 2>&1

            V1=$(awk 'FNR==3{print $3}' results/37_AAI-enveomics"$sfx"/$F1.$F2.res.out_file)
            V2=$(awk 'FNR==3{print $7}' results/37_AAI-enveomics"$sfx"/$F1.$F2.res.out_file)
            V3=$(awk 'FNR==1{print $8}' results/37_AAI-enveomics"$sfx"/$F1.$F2.res.out_file)
            V4=$(awk 'FNR==2{print $8}' results/37_AAI-enveomics"$sfx"/$F1.$F2.res.out_file)
            echo "$F1 $F2 $V1 $V2 $V3 $V4" >> results/37_AAI-enveomics"$sfx".AAI.csv
        done
    done
fi
###############################################################################
echo "completed.... ANI-enveomics --------------------------------------------"
###############################################################################
