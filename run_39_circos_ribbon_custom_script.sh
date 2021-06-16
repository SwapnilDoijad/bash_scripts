#A=1
#for F1 in $(cat list.13.txt); do
#bases=$(grep "bases: " results/08_annotation/raw_files/$F1/$F1.txt | awk '{print $2}')
#B=$(( $A + 1 ))
#echo "chr - hs$A $A 0 $bases $F1"
#A=$B
#done
(rm results/14_mugsy_alignment/all_coordinates.filtered.txt)> /dev/null 2>&1
(rm -rf results/14_mugsy_alignment/split_coordinates)> /dev/null 2>&1
(mkdir results/14_mugsy_alignment/split_coordinates/)> /dev/null 2>&1
awk '{print $1, $2, $3, $4, $5, $6}' results/14_mugsy_alignment/alignment.maf > results/14_mugsy_alignment/alignment.maf.coordinates.txt
awk 'BEGIN{i++} !NF{++i;next} {print > "filename"i".txt"}' results/14_mugsy_alignment/alignment.maf.coordinates.txt
mv filename*.txt results/14_mugsy_alignment/split_coordinates/
ls results/14_mugsy_alignment/split_coordinates/filename*.txt | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' | sed 's/filename//g' | sort -n | sed -e 's/^/filename/' > results/14_mugsy_alignment/split_coordinates/list_filenames
sed -i 1d results/14_mugsy_alignment/split_coordinates/filename1.txt

min_isolates=1
total_isolates=$( ls results/14_mugsy_alignment/fasta/*.fasta | wc -l )
#echo "total isolates = $total_isolates"
for filename in $(cat results/14_mugsy_alignment/split_coordinates/list_filenames) ; do
isolates=$(awk -F'=' 'NR==1 {print $NF}' results/14_mugsy_alignment/split_coordinates/$filename.txt)
#echo total_isolates $total_isolates isolates $isolates 
    if [ "$isolates" -gt "$min_isolates" ] && [ "$isolates" -lt "$total_isolates" ] ; then 
    echo $filename >> results/14_mugsy_alignment/split_coordinates/list_filenames_only_shared_segements  
    fi
done


for filename in $(cat results/14_mugsy_alignment/split_coordinates/list_filenames_only_shared_segements ) ; do
#for filename in $(cat results/14_mugsy_alignment/split_coordinates/list_filenames) ; do
    #thickness_tmp=$( grep score results/14_mugsy_alignment/split_coordinates/$filename.txt | awk -F'=' '{print $2}' | awk '{print $1}' )
    #thickness=$(( $thickness_tmp / 1000 ))
    #if [ $thickness -gt 80 ] ; then
    #    thickness=80
    #fi
    sed -i 1d results/14_mugsy_alignment/split_coordinates/$filename.txt
    (rm results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.txt) > /dev/null 2>&1
    for F1 in $(cat list.13.txt); do
    coordinates_1=$(grep $F1 results/14_mugsy_alignment/split_coordinates/$filename.txt | awk '{print $3, $3+$4}')
    # echo $F1 $coordinates_1
        for F2 in $(cat list.13.txt); do
            if [ "$F1" != "$F2" ] ; then
            coordinates_2=$(grep $F2 results/14_mugsy_alignment/split_coordinates/$filename.txt | awk '{print $3, $3+$4}')
            # echo $F2 $coordinates_2
            echo $F1 $coordinates_1 $F2 $coordinates_2 >> results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.txt
            fi
        done
    done 

    (results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt) > /dev/null 2>&1
    B=0
    for F1 in $(cat list.13.txt); do
    awk '$1 == "'$F1'" { print $0 }' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.txt | awk -v line=$B 'NR=='$B'' >> results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    B=$(( $B + 1 ))
    done
    #cp results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.txt results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt

    #sed -i 's/$/ thickness='$thickness'/' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/S1/hs1/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/S2/hs2/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/S3/hs3/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/S4/hs4/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/S5/hs5/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/CP017402/hs6/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/CP017403/hs7/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/CP017404/hs8/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/CP017405/hs9/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/CP034182/hs10/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/CP034101/hs11/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/CP019957/hs12/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt
    sed -i 's/NC_002929/hs13/g' results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt

    cat results/14_mugsy_alignment/split_coordinates/$filename.all_coordinates.filtered.txt | awk -F' ' '$6!=""' >> results/14_mugsy_alignment/all_coordinates.filtered.txt
done
mv  results/14_mugsy_alignment/all_coordinates.filtered.txt /home/swapnil/tools/circos-0.69/test_ribbon/all_coordinates.filtered.txt
#------------------------------------------------------------------------------
path=$(echo pwd)
#run circos
cd /home/swapnil/tools/circos-0.69
#(
./bin/circos -conf conf.ribbon.conf
#) > /dev/null 2>&1
cd $path
#------------------------------------------------------------------------------