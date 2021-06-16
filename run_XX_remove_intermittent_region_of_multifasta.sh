for F1 in $(cat list.txt); do
V1=$(cat Contamination_"$F1".txt | tail -3 | head -1 | awk -F'_' '{print $1}')
grep $V1 Contamination_"$F1".txt > $F1.2.txt

    if grep -q , "$F1.2.txt"; then
        grep ',' $F1.2.txt > $F1.comma.txt
        sed -i "\#,#d" $F1.2.txt
        awk '{print $1}' $F1.comma.txt > $F1.comma.list.txt
        rm $F1.comma.new.txt
        for C1 in $(cat $F1.comma.list.txt); do
            grep $C1 $F1.comma.txt | awk -F',' '{print $1}' $F1.comma.txt > $F1.comma.new.txt
            V2=$(grep $C1 $F1.comma.txt | awk -F',' '{print $2}' $F1.comma.txt | awk '{print $1}' )
            V3=$(grep $C1 $F1.comma.txt | awk '{print $1,$2}' $F1.comma.txt)
            echo $V3 $V2 >> $F1.comma.new.txt
        done
        cat $F1.2.txt $F1.comma.new.txt > $F1.2a.txt
        else
        cat $F1.2.txt > $F1.2a.txt
    fi
awk '{print $1, $3}' $F1.2a.txt | awk -F'.' '{print $1,$3}' OFS='\t' | sed 's/ /\t/' > $F1.3.txt
awk '{print $1}' $F1.3.txt | sort -u > $F1.list.tmp
        rm $F1.bed
        for F2 in $(cat $F1.list.tmp); do
        grep $F2 $F1.3.txt | awk '{print "'$F2'",$2,$3}' | sed 's/ /\t/g' >> $F1.bed
        done
bedtools maskfasta -fi $F1.fas -bed $F1.bed -fo $F1.filtered.fasta -mc X
sed -i -E '/>/!s/X//g' $F1.filtered.fasta    
# rm *.txt
done