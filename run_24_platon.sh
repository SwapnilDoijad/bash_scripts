###############################################################################
# 24 platon
###############################################################################
echo "started... step-24 platon -----------------------------------------------"
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/24_platon) > /dev/null 2>&1
    (mkdir results/24_platon/raw_files) > /dev/null 2>&1
###############################################################################

#for F1 in $(cat $list); do
#    (mkdir results/24_platon/raw_files/$F1) > /dev/null 2>&1 
#    /home/swapnil/tools/platon-master/bin/platon \
#    -d /media/swapnil/databases/bioinfoDBs/platon \
#    -o results/24_platon/raw_files/$F1 \
#    results/04_assembly/raw_files/$F1/$F1.contigs.fasta 
#done

#if [ ! -f results/24_platon/contigs_results.tab ]; then 
#echo "ID	Length	Coverage	# ORFs	RDS	Circular	Inc Type(s)	# Replication	# Mobilization	# OriT	# Conjugation	# AMRs	# rRNAs	# Plasmid Hits" > results/24_platon/results.tab
#fi
#for F1 in $(cat $list); do
#tail -n +2 results/24_platon/raw_files/$F1/$F1.contigs.tsv >> results/24_platon/contigs_results.tab
#done

for F1 in $(cat $list); do
contigs_sum=$( awk -F '\t' '{sum += $2} END {print sum}' results/24_platon/raw_files/$F1/$F1.contigs.tsv )
echo $F1 $contigs_sum >> results/24_platon/contigs_sum_results.tab
done
sed -i '1 i\isolate plasmid_contig_sum' results/24_platon/contigs_sum_results.tab
###############################################################################

