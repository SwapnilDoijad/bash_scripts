for F1 in (cat $list); do
hmmsearch --tblout $F1.hmm.tsv --cut_tc --notextw /home/swapnil/tools/bcgTree/data/essential.hmm $F1.faa
done
