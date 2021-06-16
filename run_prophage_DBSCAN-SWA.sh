###############################################################################
# 30 prophage dbscan
###############################################################################
echo "started... step-30 prphage dbscan -----------------------------------------"
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/30_prophage_DBSCAN-SWA) > /dev/null 2>&1
    (mkdir results/30_prophage_DBSCAN-SWA/raw_files) > /dev/null 2>&1
##############################################################################
for F1 in $(cat $list); do
echo "running prophage DBSCAN-SWA for $F1"
    dbscan-swa \
    --input results/08_annotation/raw_files/$F1/$F1.gbk \
    --output results/30_prophage_DBSCAN-SWA/raw_files/$F1 \
    --prefix $F1 \
    --add_annotation none
done
###############################################################################

