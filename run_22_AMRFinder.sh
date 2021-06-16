#!/bin/bash
###############################################################################
#22 AMRFinder
## Conda doesnt matter
## update amrfinder database
# /home/swapnil/tools/amrfinder/amrfinder.pl -u
###############################################################################
echo "started.... AMRFinder --------------------------------------------------"
###############################################################################
    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate myenv
###############################################################################
## file and directory preparation
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "want to suffix files? type "y" else press enter to continue"
    read answer
    if [ "$answer" = "y" ]; then
    s=$(echo "_$l")
    fi

    echo "Is the input is DNA (d) or protein (p)?"
    read input
    ## only Possible organisms: Campylobacter|Escherichia|Salmonella
    #echo "which organism (genus)? for e.g. Enterobacter" 
    #read organism

    (mkdir results/22_AMRFinder"$s")> /dev/null 2>&1
    (mkdir results/22_AMRFinder"$s"/results)> /dev/null 2>&1
    (mkdir results/22_AMRFinder"$s"/tmp)> /dev/null 2>&1

###############################################################################
version=$(amrfinder --version)
db_version=$(find /home/swapnil/miniconda3/envs/myenv/share/amrfinderplus/data/ -type f -exec stat \{} --printf="%y\n" \; | sort -n -r | head -n 1| awk '{print $1}') 
###############################################################################
## run amrfinder
    for F1 in $(cat $list); do
        if [ ! -f results/22_AMRFinder"$s"/results/$F1.csv ] ; then
            echo "runnning... AMRFinder $version database verssion $db_version blast for $F1"
                if [ "$input" = "p" ] ; then
                    (amrfinder -p results/08_annotation/raw_files/$F1/$F1.faa -o results/22_AMRFinder"$s"/results/$F1.csv )> /dev/null 2>&1
                    elif [ "$input" = "d" ] ; then
                    (amrfinder -n results/04_assembly/all_fasta/$F1.fasta -o results/22_AMRFinder"$s"/results/$F1.csv )> /dev/null 2>&1
                fi
                awk 'FNR>1' results/22_AMRFinder"$s"/results/$F1.csv >> results/22_AMRFinder"$s"/tmp/all.txt
            else
            echo "amrfinder already finished for $F1"
        fi
    done
###############################################################################
##calculate Abr-gene frequency
    (sed -i '/symbol/d' results/22_AMRFinder"$s"/tmp/all.txt)> /dev/null 2>&1
    heading_list=$(head -1 $list)
    heading=$(cat results/22_AMRFinder"$s"/results/$heading_list.csv | head -1)
    ex -sc "1i|$heading" -cx results/22_AMRFinder"$s"/tmp/all.txt

    cp results/22_AMRFinder"$s"/tmp/all.txt results/22_AMRFinder"$s"/tmp/all.csv # create all.csv earlier and remove this step  
    ssconvert results/22_AMRFinder"$s"/tmp/all.csv results/22_AMRFinder"$s"/tmp/all.csv.xlsx
    #------
    cp results/22_AMRFinder"$s"/tmp/all.txt results/22_AMRFinder"$s"/tmp/all.txt.1.tmp
    sed -i 's/ /_/g' results/22_AMRFinder"$s"/tmp/all.txt.1.tmp
    csvtool -t TAB -u TAB namedcol "Gene_symbol","Sequence_name","Scope","Element_type","Element_subtype","Class","Subclass" results/22_AMRFinder"$s"/tmp/all.txt.1.tmp | awk '{if(NR>1)print}' > results/22_AMRFinder"$s"/tmp/all.txt.2.tmp
    sed -i 's/\t/===/g' results/22_AMRFinder"$s"/tmp/all.txt.2.tmp
    cat results/22_AMRFinder"$s"/tmp/all.txt.2.tmp | sed 's/===/=/g' | sort -t '=' -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 | sed 's/=/===/g' | uniq > results/22_AMRFinder"$s"/tmp/all.txt.3.tmp 

    (rm results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv.tmp)> /dev/null 2>&1
    for F2 in $(cat results/22_AMRFinder"$s"/tmp/all.txt.3.tmp); do
    V1=$(grep -c "$F2" results/22_AMRFinder"$s"/tmp/all.txt.2.tmp)
    #echo $V1 $F2 >> results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv.tmp
    echo $V1 $F2 >> results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv
    done

    #cat results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv.tmp | sort -k 1,1rn results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv.tmp > results/22_AMRFinder"$s"/Abr_gene_frequency.csv

    ex -sc '1i|frequency Abr-gene description Scope Element_type Element_subtype Class Subclass' -cx results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv

    sed -i 's/===/ /g' results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv
    ssconvert results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv results/22_AMRFinder"$s"/Abr_gene_frequency.csv.xlsx
    #rm results/22_AMRFinder"$s"/tmp/Abr_gene_frequency.csv.tmp

###############################################################################
## Creating 1-0 matrix of ABR results
    echo "running matrix for ABR-DB"
    awk -F'===' '{print $1}' results/22_AMRFinder"$s"/tmp/all.txt.3.tmp > results/22_AMRFinder"$s"/tmp/all.txt.genes.tmp

    for F1 in $(cat $list);do
    csvtool -t TAB -u TAB namedcol "Gene symbol" results/22_AMRFinder"$s"/results/$F1.csv | awk '{if(NR>1)print}' | sed 's/ /_/g' > results/22_AMRFinder"$s"/tmp/$F1.gene-list-to-count.tmp1
        (rm results/22_AMRFinder"$s"/tmp/$F1.Abr-gene-count.tmp1) > /dev/null 2>&1 
        for V1 in $(cat results/22_AMRFinder"$s"/tmp/all.txt.genes.tmp); do
        V2=$(awk '{count[$1]++} END {print count["'$V1'"]}' results/22_AMRFinder"$s"/tmp/$F1.gene-list-to-count.tmp1)
        echo $V2 >> results/22_AMRFinder"$s"/tmp/$F1.Abr-gene-count.tmp1
        done
    awk '{for (i=1; i<= NF; i++) {if($i > 1) { $i=1; } } print }' results/22_AMRFinder"$s"/tmp/$F1.Abr-gene-count.tmp1 > results/22_AMRFinder"$s"/tmp/$F1.Abr-gene-count.tmp2
    ex -sc '1i|'$F1'' -cx results/22_AMRFinder"$s"/tmp/$F1.Abr-gene-count.tmp1
    ex -sc '1i|'$F1'' -cx results/22_AMRFinder"$s"/tmp/$F1.Abr-gene-count.tmp2
    sed -i -e 's/^$/0/' results/22_AMRFinder"$s"/tmp/$F1.Abr-gene-count.tmp1
    sed -i -e 's/^$/0/' results/22_AMRFinder"$s"/tmp/$F1.Abr-gene-count.tmp2
    done
    
    cp results/22_AMRFinder"$s"/tmp/all.txt.genes.tmp results/22_AMRFinder"$s"/tmp/all.txt.4.2.tmp
    sed -i 's/$/:c/' results/22_AMRFinder"$s"/tmp/all.txt.4.2.tmp
    sed -i '1i\===\' results/22_AMRFinder"$s"/tmp/all.txt.4.2.tmp

    ls -1 results/22_AMRFinder"$s"/tmp/*.Abr-gene-count.tmp1 | split -l 500 -d - lists1
    for list1 in lists1* ; do
        paste $(cat $list1) > results/22_AMRFinder"$s"/tmp/merge.Abr-gene-count.tmp1.${list1##lists1}; 
    done
    paste results/22_AMRFinder"$s"/tmp/merge.Abr-gene-count.tmp1.* > results/22_AMRFinder"$s"/tmp/all.Abr-gene-count.1.tab
    rm results/22_AMRFinder"$s"/tmp/*.Abr-gene-count.tmp1.*
    rm lists1*

    ls -1 results/22_AMRFinder"$s"/tmp/*.Abr-gene-count.tmp2  | split -l 500 -d - lists2
    for list2 in lists* ; do
        paste $(cat $list2) > results/22_AMRFinder"$s"/tmp/merge.Abr-gene-count.tmp2.${list2##lists2}; 
    done
    paste results/22_AMRFinder"$s"/tmp/merge.Abr-gene-count.tmp2.* > results/22_AMRFinder"$s"/tmp/all.Abr-gene-count.2.tab
    rm results/22_AMRFinder"$s"/tmp/*.Abr-gene-count.tmp2.*
    rm lists2*

    paste results/22_AMRFinder"$s"/tmp/all.txt.4.2.tmp results/22_AMRFinder"$s"/tmp/all.Abr-gene-count.1.tab > results/22_AMRFinder"$s"/tmp/matrix.csv.tmp1
    paste results/22_AMRFinder"$s"/tmp/all.txt.4.2.tmp results/22_AMRFinder"$s"/tmp/all.Abr-gene-count.2.tab > results/22_AMRFinder"$s"/tmp/matrix_1-0.csv.tmp2

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
    }' results/22_AMRFinder"$s"/tmp/matrix.csv.tmp1 > results/22_AMRFinder"$s"/matrix.csv

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
    }' results/22_AMRFinder"$s"/tmp/matrix_1-0.csv.tmp2 > results/22_AMRFinder"$s"/matrix_1-0.csv
###############################################################################
## Create final files
    sed -i 's/===//g' results/22_AMRFinder"$s"/matrix.csv
    sed -i 's/===//g' results/22_AMRFinder"$s"/matrix_1-0.csv
    sed -i 's/ /,/g' results/22_AMRFinder"$s"/matrix_1-0.csv 
    echo "finished matrix for ABR-DB"
###############################################################################
echo "Finished.... AMRFinder  ------------------------------------------------"
###############################################################################
    conda deactivate
###############################################################################