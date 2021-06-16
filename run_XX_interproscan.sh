## ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/lookup_service/ 
## https://interproscan5-docs.readthedocs.io/en/latest/LocalLookupService.html

(mkdir results/08_annotation/interproscan_files) > /dev/null 2>&1
for F1 in $(cat list.all.txt); do
/home/swapnil/tools/interproscan-5.19-58.0/interproscan.sh \
-i results/08_annotation/raw_files/$F1/$F1.faa \
-iprlookup \
-pa \
-f tsv \
-o results/08_annotation/interproscan_files/$F1.tab
done