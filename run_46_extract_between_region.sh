## either provide location or provide reference_fasta
###############################################################################
## extract upstream downstrema region
###############################################################################
## file preparation
    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate myenv

    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "provide upstream reference fasta (for e.g. data/genes/arn_region.fasta)"
    echo "-------------------------------------------------------------------------"
    ls data/genes/*.fasta | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read region_up
    echo "region_up $region_up"

    echo "provide downstream reference fasta (for e.g. data/genes/arn_region.fasta)"
    echo "-------------------------------------------------------------------------"
    ls data/genes/*.fasta | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read region_down
    echo "region_down $region_down"

    region_name1=$(echo "$region_up" | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' )
    region_name2=$(echo "$region_down" | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' )
    region_name=$(echo $region_name1 $region_name2 | sed 's/ /_/g')
    echo $region_name

    echo "need to annotate the region? press y else press ENTER"
    read annot

    (mkdir results/46_extract_between_"$region_name")> /dev/null 2>&1
    path1=results/46_extract_between_"$region_name"
    (mkdir $path1)> /dev/null 2>&1
    (mkdir $path1/tmp)> /dev/null 2>&1
    (mkdir $path1/fasta)> /dev/null 2>&1
    (mkdir $path1/fasta_extracted)> /dev/null 2>&1
    (mkdir $path1/fasta_extracted/discarded)> /dev/null 2>&1
    (mkdir $path1/annotations)> /dev/null 2>&1
    (mkdir $path1/annotations/raw_files)> /dev/null 2>&1
    (mkdir $path1/annotations/gbk)> /dev/null 2>&1
###############################################################################
for F1 in $(cat $list); do
    cp results/04_assembly/all_fasta/$F1.fasta $path1/fasta/$F1.fasta

    ## format fasta 
    sed -i 's/>.*/NNNNNNNNNN/g' $path1/fasta/$F1.fasta
    sed -i '0,/NNNNNNNNNN/s///' $path1/fasta/$F1.fasta
    sed -i "1i "'>'$F1"" $path1/fasta/$F1.fasta    
    sed -i '/^[[:space:]]*$/d' $path1/fasta/$F1.fasta 
    (python /home/swapnil/tools/read-cleaning-format-conversion/KSU_bioinfo_lab/fasta-o-matic/fasta_o_matic.py -o $path1/fasta/ -f $path1/fasta/$F1.fasta)> /dev/null 2>&1
    (mv $path1/fasta/"$F1"_wrap.fasta $path1/fasta/$F1.fasta)> /dev/null 2>&1
    
    ## find co-ordinates of up gene ($region_name1)
    blastn -subject $path1/fasta/$F1.fasta -query $region_up -out $path1/tmp/$F1.custom-gene-blast.up.tmp -max_target_seqs 1 -max_hsps 1 -evalue 1e-10 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
    if [ -s "$path1/tmp/$F1.custom-gene-blast.up.tmp" ]; then
        up_start=$(awk 'NR==1 {print $3}' $path1/tmp/$F1.custom-gene-blast.up.tmp)
        up_end=$(awk 'NR==1 {print $4}' $path1/tmp/$F1.custom-gene-blast.up.tmp)
        else
        echo "$F1 $region_name1 not found"
    fi
    ## find co-ordinates of down gene ($region_name2)
    blastn -subject $path1/fasta/$F1.fasta -query $region_down -out $path1/tmp/$F1.custom-gene-blast.down.tmp -max_target_seqs 1 -max_hsps 1 -evalue 1e-10 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
    if [ -s "$path1/tmp/$F1.custom-gene-blast.down.tmp" ]; then
        down_start=$(awk 'NR==1 {print $3}' $path1/tmp/$F1.custom-gene-blast.down.tmp)
        down_end=$(awk 'NR==1 {print $4}' $path1/tmp/$F1.custom-gene-blast.down.tmp)
        else
        echo "$F1 $region_name2 not found"
    fi

    coordinates=$(echo $up_end $up_start $down_end $down_start | tr " " "\n" | sort -u )
    #echo $F1 $coordinates | tr " " "\n"
    final_up=$(echo $coordinates | tr " " "\n" | head -1)
    final_up_100=$(( $final_up - 100 ))
    final_down=$(echo $coordinates | tr " " "\n" | tail -1)
    final_down_100=$(( $final_down + 100 ))
    region_length=$(( $final_down - $final_up ))
    #echo $final_up $final_down
    echo "$F1: $region_name is $region_length bp long"
    echo "$F1: $region_name is $region_length bp long" >> $path1/results.csv

    (rm $path1/fasta/$F1.fasta.fai)> /dev/null 2>&1
    echo -e "$F1\t$final_up_100\t$final_down_100\t"$F1"_region" > $path1/tmp/$F1.bed
    (bedtools getfasta -fi $path1/fasta/$F1.fasta -bed $path1/tmp/$F1.bed -name > $path1/fasta_extracted/"$F1"_region.fasta)> /dev/null 2>&1

    #find if reverse complement and align as per reference
    if [ "$up_start" -gt "$down_start" ] ; then
    sed '1d' $path1/fasta_extracted/"$F1"_region.fasta | perl -pe 'chomp;tr/ACGTNacgtn/TGCANtgcan/;$_=reverse."\n"' | sed '1i >'$F1'' > $path1/tmp/"$F1"_region.fasta
    mv $path1/tmp/"$F1"_region.fasta $path1/fasta_extracted/"$F1"_region.fasta
    fi

    sed -i '1d' $path1/fasta_extracted/"$F1"_region.fasta 
    sed -i '1i >'$F1"_$region_name"'' $path1/fasta_extracted/"$F1"_region.fasta

    ## annotate the region
    if [ "$annot" == "n" ]; then
    :
    else
    (prokka --quiet --outdir $path1/annotations/raw_files/$F1 --force --prefix $F1 --locustag $F1 --strain $F1 --rnammer $path1/fasta_extracted/"$F1"_region.fasta)> /dev/null 2>&1
    cp $path1/annotations/raw_files/$F1/$F1.gbk $path1/annotations/gbk/$F1.gbk
    fi

done
conda deactivate

###############################################################################

###############################################################################