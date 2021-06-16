#!/bin/bash
###############################################################################
#11 ANI
###############################################################################
echo "started... step-11 ANI --------------------------------------------------"
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")



(mkdir results)> /dev/null 2>&1
(mkdir results/11_ANI)> /dev/null 2>&1
echo "query subject ANI alinged(bp) 1020-fragments-aligned %-genome-covered" > results/11_ANI/ANI.csv
###############################################################################
for F1 in $(cat $list); do
echo "running.... step-11 ANI for..." $F1 

    for F2 in $(cat $list); do
    echo "running $F1 vs $F2"
        (perl /home/swapnil/pipeline/tools/ANI_and_Total-aligned_plasmid_200.pl \
        --fd /home/swapnil/tools/blast-2.2.9-amd64-linux/bin/formatdb \
        --bl /home/swapnil/tools/blast-2.2.9-amd64-linux/bin/blastall \
        --qr results/04_assembly/all_fasta/$F1.fasta \
        --sb results/04_assembly/all_fasta/$F2.fasta \
        --od results/11_ANI/output \
        >> results/11_ANI/ANI.csv)> /dev/null 2>&1
    done

echo "finished... step-11 ANI for..." $F1
done

rm error.log formatdb.log 

sed -i 's/\/home\/swapnil\/test\/results\/04_assembly\/all_fasta\///g' results/11_ANI/ANI.csv
sed -i 's/\.fasta//g' results/11_ANI/ANI.csv

###############################################################################
## matrix-creation step

for F1 in $(cat $list);do

	for F2 in $(cat $list); do
        V1=$(awk '$1 == "'$F1'" && $2 == "'$F2'" {print 100-$3}' results/11_ANI/ANI.csv)
        if [ -z $V1 ] ; then
            echo 100 >> results/11_ANI/output/$F1.ANI.tmp
            else
            echo $V1 >> results/11_ANI/output/$F1.ANI.tmp
        fi
    done
    sed -i '1 i\'$F1'' results/11_ANI/output/$F1.ANI.tmp
done

touch results/11_ANI/output/all.ANI.tmp
for F3 in $(cat $list); do
paste results/11_ANI/output/all.ANI.tmp results/11_ANI/output/$F3.ANI.tmp > results/11_ANI/output/all.ANI.tmp.tmp
mv results/11_ANI/output/all.ANI.tmp.tmp results/11_ANI/output/all.ANI.tmp 
done
awk 'NR==1{print}' results/11_ANI/output/all.ANI.tmp | /home/swapnil/pipeline/tools/transpose.sh > results/11_ANI/output/isolate-name.txt
sed -i '1s/^/\n/' results/11_ANI/output/isolate-name.txt
paste results/11_ANI/output/isolate-name.txt results/11_ANI/output/all.ANI.tmp > results/11_ANI/all.ANI.tab
Rscript /home/swapnil/pipeline/tools/plot_and_tree.fastANI-distance.r
mv Rplots.pdf results/11_ANI/
mv tree.nwk results/11_ANI/

#------------------------------------------------------------------------------
# get a representative of a clutster at 95% (-h 5)
echo "running ggrasp.R"
/home/swapnil/tools/GGRaSP-master/ggrasp.R -i results/11_ANI/all.ANI.tab -d 100 -h 5 -o results/11_ANI/representative_of_a_clusters_at_95_ANI --plottree 
rm Rplots.pdf
###############################################################################
echo "completed... step-11 ANI ------------------------------------------------"
###############################################################################
