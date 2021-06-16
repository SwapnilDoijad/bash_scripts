## either provide location or provide reference_fasta

###############################################################################
## extract upstream downstrema region
###############################################################################
echo "provide list file (for e.g. all)"
echo "-------------------------------------------------------------------------"
ls list.*.txt | sed 's/ /\n/g'
echo "-------------------------------------------------------------------------"
read l
list=$(echo "list.$l.txt")

echo "provide reference fasta (for e.g. data/genes/arn_region.fasta)"
echo "-------------------------------------------------------------------------"
ls data/genes/*.fasta | sed 's/ /\n/g'
echo "-------------------------------------------------------------------------"
read region

region_name=$(echo $region | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}')

echo how much bp upstream and downstream?
read UpDown

echo "need to annotate the region? press y else press ENTER"
read annot

path1=results/46_extract_U_D_"$region_name"
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
    echo "running $F1"
    ## check if closed version of fasta if available and provide path for the file
    V1=$(find -name $F1.fasta | grep 'DnaA_arranged')
    #if [ ! -z "$V1" ]; then
    #    cp /home/swapnil/p_Entb_Germany/results_closed/08_annotation/raw_files/$F1/$F1.fna  $path1/fasta/$F1.fasta
    #    else
    #    cp /home/swapnil/p_Entb_Germany/results/08_annotation/raw_files/$F1/$F1.fna $path1/fasta/$F1.fasta
    #fi
    cp results/04_assembly/all_fasta/$F1.fasta $path1/fasta/$F1.fasta

    ## format fasta 
    sed -i 's/>.*/NNNNNNNNNN/g' $path1/fasta/$F1.fasta
    sed -i '0,/NNNNNNNNNN/s///' $path1/fasta/$F1.fasta
    sed -i "1i "'>'$F1"" $path1/fasta/$F1.fasta    
    sed -i '/^[[:space:]]*$/d' $path1/fasta/$F1.fasta 
    (python /home/swapnil/tools/read-cleaning-format-conversion/KSU_bioinfo_lab/fasta-o-matic/fasta_o_matic.py -o $path1/fasta/ -f $path1/fasta/$F1.fasta)> /dev/null 2>&1
    (mv $path1/fasta/"$F1"_wrap.fasta $path1/fasta/$F1.fasta)> /dev/null 2>&1
    
    blastn -subject $path1/fasta/$F1.fasta -query $region -out $path1/tmp/$F1.custom-gene-blast.tmp -max_target_seqs 1 -max_hsps 1 -evalue 1e-10 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
    if [ -s "$path1/tmp/$F1.custom-gene-blast.tmp" ]; then
        start=$(awk 'NR==1 {print $3}' $path1/tmp/$F1.custom-gene-blast.tmp)
        end=$(awk 'NR==1 {print $4}' $path1/tmp/$F1.custom-gene-blast.tmp)
        if [ "$start" -gt "$end" ]; then
            Forw=$(( $end - $UpDown ))
            Revr=$(( $start + $UpDown))
            (rm $path1/fasta/$F1.fasta.fai)> /dev/null 2>&1
            echo -e "$F1\t$Forw\t$Revr\t"$F1"_region" > $path1/tmp/$F1.bed
            (bedtools getfasta -fi $path1/fasta/$F1.fasta -bed $path1/tmp/$F1.bed -name > $path1/fasta_extracted/"$F1"_region.RC.fasta)> /dev/null 2>&1
            sed '1d' $path1/fasta_extracted/"$F1"_region.RC.fasta | perl -pe 'chomp;tr/ACGTNacgtn/TGCANtgcan/;$_=reverse."\n"' | sed '1i >'$F1'' > $path1/fasta_extracted/"$F1"_region.fasta
            mv $path1/fasta_extracted/"$F1"_region.RC.fasta $path1/fasta_extracted/discarded/
            else
            Forw=$(( $start - $UpDown ))
            Revr=$(( $end + $UpDown ))
            (rm $path1/fasta/$F1.fasta.fai)> /dev/null 2>&1
            echo -e "$F1\t$Forw\t$Revr\t"$F1"_region" > $path1/tmp/$F1.bed
            (bedtools getfasta -fi $path1/fasta/$F1.fasta -bed $path1/tmp/$F1.bed -name > $path1/fasta_extracted/"$F1"_region.fasta)> /dev/null 2>&1
        fi

        ## annotate the region
        if [ "$annot" == "y" ]; then
        prokka --quiet --outdir $path1/annotations/raw_files/$F1 --force --prefix $F1 --locustag $F1 --strain $F1 --rnammer $path1/fasta_extracted/"$F1"_region.fasta
        cp $path1/annotations/raw_files/$F1/$F1.gbk $path1/annotations/gbk/$F1.gbk
        fi

        echo "$F1" >> $path1/region_present_in_strains.list
        else
        echo "$F1" >> $path1/region_absent_in_strains.list
    fi
done
###############################################################################

###############################################################################