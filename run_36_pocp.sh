#!/bin/bash
###############################################################################
#36 POCP (percentage of conserved protein)

###############################################################################

echo "started.... POCP -------------------------------------------------------"

###############################################################################
if [ ! -d results/36_POCP ]; then
mkdir results/36_POCP
fi

if [ ! -d results/36_POCP/tmp ]; then
mkdir results/36_POCP/tmp
fi

rm results/36_POCP/number-of-conserved-proteins.tab > /dev/null 2>&1

for F1 in $(cat list.txt);do
V1=$(grep -c ">" results/08_annotation/all_faa/$F1.faa)

if [ ! -d results/36_POCP/tmp/$F1 ]; then
mkdir results/36_POCP/tmp/$F1
fi

for F2 in $(cat list.txt);do
V2=$(grep -c ">" results/08_annotation/all_faa/$F2.faa)

blastp -subject results/08_annotation/all_faa/$F1.faa -query results/08_annotation/all_faa/$F2.faa -out results/36_POCP/tmp/$F1/$F1.$F2.result.tmp -max_hsps 1 -max_target_seqs 1 -evalue 1e-5 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"

awk '$14 >= 40 && $15 >= 50 {print $0}' results/36_POCP/tmp/$F1/$F1.$F2.result.tmp > results/36_POCP/tmp/$F1/$F1.$F2.result.tab ;

V3=$(wc -l results/36_POCP/tmp/$F1/$F1.$F2.result.tab | awk '{print $1}' )

#number-of-proteins-in-subject number-of-proteins-in-query and number-of-shared-proteins(>40% protein identity and >50% coverage)
echo "$V1 $V2 $V3" > results/36_POCP/tmp/$F1/$F1.$F2.number-of-shared-proteins.tab

done
done

for F1 in $(cat list.txt);do
V1=$(grep -c ">" results/08_annotation/all_faa/$F1.faa)
for F2 in $(cat list.txt);do
V2=$(grep -c ">" results/08_annotation/all_faa/$F2.faa)

V4=$(awk '{print $3}' results/36_POCP/tmp/$F1/$F1.$F2.number-of-shared-proteins.tab)
V5=$(awk '{print $3}' results/36_POCP/tmp/$F2/$F2.$F1.number-of-shared-proteins.tab)
V6=$(($V4+$V5))
V7=$(($V1+$V2))

V8=$(echo "scale=2; $V6*100/$V7" | bc -l)

echo "the percentage of conserved protein (POCP) between $F1 and $F2 = $V8" >> results/36_POCP/number-of-conserved-proteins.tab

done
done


###############################################################################

echo "completed.... POCP -----------------------------------------------------"

###############################################################################
