#!/bin/bash
###############################################################################
#33 minhash
    ## the mash distance can also be preapred by 'mash dist mash_sketch.msh mash_sketch.msh| square_mash > mash.tsv'
    ## Need to incorporate this
###############################################################################
echo "started.... step-33 minhash --------------------------------------------"
###############################################################################
## file and directory preparatoins
    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")
    #------------------------------------------------------------------------------
    (mkdir results/33_minhash_"$l")> /dev/null 2>&1
    (mkdir results/33_minhash_"$l"/fasta)> /dev/null 2>&1
    (mkdir results/33_minhash_"$l"/tmp)> /dev/null 2>&1
    (mkdir results/33_minhash_"$l"/msh)> /dev/null 2>&1
    (mkdir results/33_minhash_"$l"/msh/distance)> /dev/null 2>&1
    (mkdir results/33_minhash_"$l"/msh/distance/tmp)> /dev/null 2>&1
###############################################################################
## Mash with r package
    echo "copying fasta files"
    for F1 in $(cat $list); do
        cp results/04_assembly/all_fasta/"$F1".fasta results/33_minhash_"$l"/fasta/
    done
    echo "sketching"
    (mash sketch -o results/33_minhash_"$l"/msh/mash_sketch.msh results/33_minhash_"$l"/fasta/*.fasta)> /dev/null 2>&1
    rm -rf results/33_minhash_"$l"/fasta
    echo "mashing"
    mash dist results/33_minhash_"$l"/msh/mash_sketch.msh results/33_minhash_"$l"/msh/mash_sketch.msh | square_mash > results/33_minhash_"$l"/tmp/all.distance.for_R.tab

    R1_rows=$(wc -l results/33_minhash_"$l"/tmp/all.distance.for_R.tab | awk '{print $1}' )
    R2_columns=$(head -n1 results/33_minhash_"$l"/tmp/all.distance.for_R.tab |  sed 's/\t/\n/g' | wc -l)
    echo "rows:$R1_rows columns:$R2_columns"
    #------------------------------------------------------------------------------
    # for_R
    if [ ! "$R1_rows" == "$R2_columns" ] ; then
    echo "rows and cloumns numbers are not OK for R, could not run plot_and_tree.mash-distance.r"
    else
    echo "running R"
    sed -i "s/minhash/minhash_"$l"/g" /home/swapnil/pipeline/tools/plot_and_tree.mash-distance.r
    Rscript /home/swapnil/pipeline/tools/plot_and_tree.mash-distance.r
    sed -i "s/minhash_"$l"/minhash/g" /home/swapnil/pipeline/tools/plot_and_tree.mash-distance.r
    (mv Rplots.pdf results/33_minhash_"$l"/)> /dev/null 2>&1
    (mv tree.nwk results/33_minhash_"$l"/)> /dev/null 2>&1
    fi

###############################################################################
## Mash with GGRaSP package 
    echo "running GGRaSP package"
    for F1 in $(cat $list); do
    F2=$(echo $F1 | awk -F'.' '{print $1}')
    awk -F'\t' -v col="$F2" 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print (1-$c)*100}' results/33_minhash_"$l"/tmp/all.distance.for_R.tab > results/33_minhash_"$l"/msh/distance/tmp/$F1.ggrasp.tmp
    done
    ## get first-name column
    awk -F'\t' '{print $1}' results/33_minhash_"$l"/tmp/all.distance.for_R.tab | sed '1d' > results/33_minhash_"$l"/tmp/list.ggrasp

    cat results/33_minhash_"$l"/tmp/list.ggrasp | sed '1i\\' > results/33_minhash_"$l"/msh/distance/all.ggrasp.tab
    for F1 in $(cat results/33_minhash_"$l"/tmp/list.ggrasp); do
    paste results/33_minhash_"$l"/msh/distance/all.ggrasp.tab results/33_minhash_"$l"/msh/distance/tmp/$F1.ggrasp.tmp > results/33_minhash_"$l"/msh/distance/all.ggrasp.tab.tmp
    mv results/33_minhash_"$l"/msh/distance/all.ggrasp.tab.tmp results/33_minhash_"$l"/msh/distance/all.ggrasp.tab
    done

    mv results/33_minhash_"$l"/msh/distance/all.ggrasp.tab results/33_minhash_"$l"/tmp/all.distance.for_GGRaSP.tab
    #------------------------------------------------------------------------------
    G1_rows=$(wc -l results/33_minhash_"$l"/tmp/all.distance.for_GGRaSP.tab | awk '{print $1}' )
    G2_columns=$(head -n1 results/33_minhash_"$l"/tmp/all.distance.for_GGRaSP.tab |  sed 's/\t/\n/g' | wc -l)
    echo "$G1_rows $G2_columns "
    #------------------------------------------------------------------------------
    # for_GGRaSP get a representative of a clutster at 95% (-h 5)
    if [ ! "$G1_rows" == "$G2_columns" ] ; then
    echo "rows and cloumns numbers are not OK for GGRaSP, could not run all.distance.for_GGRaSP.tab"
    else
    echo "running ggrasp.R"
    (/home/swapnil/tools/GGRaSP-master/ggrasp.R -i results/33_minhash_"$l"/tmp/all.distance.for_GGRaSP.tab -d 100 -h 5.0000 -o results/33_minhash_"$l"/representative_of_a_clusters_at_95.0000_250000nucl-5Mb-genome_ANI )> /dev/null 2>&1
    (/home/swapnil/tools/GGRaSP-master/ggrasp.R -i results/33_minhash_"$l"/tmp/all.distance.for_GGRaSP.tab -d 100 -h 0.0100 -o results/33_minhash_"$l"/representative_of_a_clusters_at_99.9900_000500nucl-5Mb-genome_ANI )> /dev/null 2>&1
    (/home/swapnil/tools/GGRaSP-master/ggrasp.R -i results/33_minhash_"$l"/tmp/all.distance.for_GGRaSP.tab -d 100 -h 0.0050 -o results/33_minhash_"$l"/representative_of_a_clusters_at_99.9950_000250nucl-5Mb-genome_ANI )> /dev/null 2>&1
    (/home/swapnil/tools/GGRaSP-master/ggrasp.R -i results/33_minhash_"$l"/tmp/all.distance.for_GGRaSP.tab -d 100 -h 0.0010 -o results/33_minhash_"$l"/representative_of_a_clusters_at_99.9990_000050nucl-5Mb-genome_ANI )> /dev/null 2>&1
    (/home/swapnil/tools/GGRaSP-master/ggrasp.R -i results/33_minhash_"$l"/tmp/all.distance.for_GGRaSP.tab -d 100 -h 0.0001 -o results/33_minhash_"$l"/representative_of_a_clusters_at_99.9999_000005nucl-5Mb-genome_ANI )> /dev/null 2>&1
    fi
    
    #----------------------------------------------------------------------------
    ## for 
        ls results/33_minhash_"$l"//*.medoids.txt | awk -F'/' '{print $NF}' > results/33_minhash_"$l"/tmp/list.ggrasp.tmp
            for F3 in $(cat results/33_minhash_"$l"/tmp/list.ggrasp.tmp ); do
                echo "running mash and R for $F3"
                name=$( echo $F3 | awk -F'/' '{print $NF}' )
                mkdir results/33_minhash_"$l"/fasta_$F3
                for files in $(cat results/33_minhash_"$l"/$F3 ); do
                    cp results/04_assembly/all_fasta/"$files".fasta results/33_minhash_"$l"/fasta_$F3/
                done
                #echo "sketching"
                (mash sketch -o results/33_minhash_"$l"/msh/mash_sketch.$name.msh results/33_minhash_"$l"/fasta_$F3/*.fasta)> /dev/null 2>&1
                rm -rf results/33_minhash_"$l"/fasta_$F3
                #echo "mashing"
                mash dist results/33_minhash_"$l"/msh/mash_sketch.$name.msh results/33_minhash_"$l"/msh/mash_sketch.$name.msh | square_mash > results/33_minhash_"$l"/tmp/all.distance.for_R.$name.tab

                #------------------------------------------------------------------------------
                # for_R
                sed -i "s/minhash/minhash_"$l"/g" /home/swapnil/pipeline/tools/plot_and_tree.mash-distance.r
                sed -i "s/all.distance.for_R.tab/all.distance.for_R.$name.tab/g" /home/swapnil/pipeline/tools/plot_and_tree.mash-distance.r
                Rscript /home/swapnil/pipeline/tools/plot_and_tree.mash-distance.r
                sed -i "s/minhash_"$l"/minhash/g" /home/swapnil/pipeline/tools/plot_and_tree.mash-distance.r
                sed -i "s/all.distance.for_R.$name.tab/all.distance.for_R.tab/g" /home/swapnil/pipeline/tools/plot_and_tree.mash-distance.r
                (mv Rplots.pdf results/33_minhash_"$l"/$name.Rplots.pdf)> /dev/null 2>&1
                (mv tree.nwk results/33_minhash_"$l"/$name.tree.nwk)> /dev/null 2>&1
            done
exit
###############################################################################
## thinspace
## did not work as exepcted so not added




###############################################################################
# screening for the closest strain (time consuming process)
    echo "Do you want to find the closest genome in NCBI? yes or PRESS enter to skip the screening"
    read answer

    if [ "$answer" == "yes" ] ; then
    for F1 in $(cat $list);do
    echo screening "$F1"
    (mash screen -w -p 8 /home/swapnil/pipeline/tools/databases/mash/refseq.genomes.k21s1000.msh results/33_minhash_"$l"/fasta/"$F1".fasta > results/33_minhash_"$l"/tmp/"$F1".screen.tmp)> /dev/null 2>&1
    sort -gr results/33_minhash_"$l"/tmp/"$F1".screen.tmp | head > results/33_minhash_"$l"/results/"$F1".screen.tab
    sed -i 's/^/'"$F1"'	/' results/33_minhash_"$l"/results/"$F1".screen.tab
    (mash screen -w -p 8 /home/swapnil/pipeline/tools/databases/mash/refseq.plasmid.k21s1000.msh results/33_minhash_"$l"/fasta/"$F1".fasta > results/33_minhash_"$l"/tmp/"$F1".screen.plasmid.tmp)> /dev/null 2>&1
    sort -gr results/33_minhash_"$l"/tmp/"$F1".screen.plasmid.tmp | head > results/33_minhash_"$l"/results/"$F1".screen.plasmid.tab
    sed -i 's/^/'"$F1"'	/' results/33_minhash_"$l"/results/"$F1".screen.plasmid.tab
    done
    #------
    cat results/33_minhash_"$l"/results/*.screen.tab > results/33_minhash_"$l"/results/all.screen.tab

    if grep -Fxq  "Isolate-id	identity	shared-hashes	median-multiplicity	p-value	query-ID	query-comment" results/33_minhash_"$l"/results/all.screen.tab; then
    :
    else
    ex -sc '1i|Isolate-id	identity	shared-hashes	median-multiplicity	p-value	query-ID	query-comment' -cx results/33_minhash_"$l"/results/all.screen.tab
    fi

    unoconv -i FilterOptions=09,,system,1 -f xls -o results/33_minhash_"$l"/all.screen.xls results/33_minhash_"$l"/results/all.screen.tab 

    #------

    cat results/33_minhash_"$l"/results/*.screen.plasmid.tab > results/33_minhash_"$l"/results/all.screen.plasmid.tab

    if grep -Fxq  "Isolate-id	identity	shared-hashes	median-multiplicity	p-value	query-ID	query-comment" results/33_minhash_"$l"/results/all.screen.plasmid.tab; then
    :
    else
    ex -sc '1i|Isolate-id	identity	shared-hashes	median-multiplicity	p-value	query-ID	query-comment' -cx results/33_minhash_"$l"/results/all.screen.plasmid.tab
    fi

    unoconv -i FilterOptions=09,,system,1 -f xls -o results/33_minhash_"$l"/all.screen.plasmid.xls results/33_minhash_"$l"/results/all.screen.plasmid.tab 

    fi
###############################################################################
echo "completed.... step-33 minhash ------------------------------------------"
###############################################################################
