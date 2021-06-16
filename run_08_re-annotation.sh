#!/bin/bash
###############################################################################

# 08 re-annotation

###############################################################################

echo "started... step-08 re-annotation ---------------------------------------"

###############################################################################
echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

(mkdir results/08_annotation/00_re-annotated) > /dev/null 2>&1
DirPth=$(echo "results/08_annotation/00_re-annotated")

for F1 in $(cat $list); do
    (mkdir $DirPth/$F1) > /dev/null 2>&1
    cp results/08_annotation/raw_files/$F1/$F1.gbk $DirPth/$F1/$F1.re-annotated.gbk

    ##-----------------------------------------------------------------------------
    ## transposases
    echo "running transposases annotation for $F1"
    csvcut -c Start,End,locus_tag,product results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/,/\t/g' > $DirPth/$F1/$F1.transposon.list.tmp
    grep -i 'transposon' $DirPth/$F1/$F1.transposon.list.tmp | awk '{print $3}' > $DirPth/$F1/$F1.transposon.list
    n=$(wc -l $DirPth/$F1/$F1.transposon.list | awk '{print $1}')
    if [ $n -gt 1 ] ; then
        for trans in $(cat $DirPth/$F1/$F1.transposon.list);do
            start=$(grep $trans $DirPth/$F1/$F1.transposon.list.tmp | awk 'NR==1 {print $1}' )
            end=$(grep $trans $DirPth/$F1/$F1.transposon.list.tmp | awk 'NR==1 {print $2}' )
            n1=$(grep -n "ORIGIN" $DirPth/$F1/$F1.re-annotated.gbk | awk -F':' '{print $1}')
            n2=$(( $n1 - 1))
            V1=$(echo '&&&&&misc_feature&&&&'$start'..'$end'==&&&&&&&&&&&&&&&&&&&&&/ugene_group='MGE'==&&&&&&&&&&&&&&&&&&&&&/ugene_name='MGE'==&&&&&&&&&&&&&&&&&&&&&/gene='transposases'==&&&&&&&&&&&&&&&&&&&&&/note='transposases'' | sed "s/ //g")
            sed -i -e ''$n2'a\'$V1'' $DirPth/$F1/$F1.re-annotated.gbk
        done
    fi

    ##-----------------------------------------------------------------------------
    ## integrases
    echo "running integrases annotation for $F1"
    csvcut -c Start,End,locus_tag,product results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/,/\t/g' > $DirPth/$F1/$F1.integrase.list.tmp
    grep -i 'integrase' $DirPth/$F1/$F1.transposon.list.tmp | awk '{print $3}' > $DirPth/$F1/$F1.integrase.list
    n=$(wc -l $DirPth/$F1/$F1.integrase.list | awk '{print $1}')
    if [ $n -gt 1 ] ; then
        for inte in $(cat $DirPth/$F1/$F1.integrase.list); do
            start=$(grep $inte $DirPth/$F1/$F1.integrase.list.tmp | awk 'NR==1 {print $1}' )
            end=$(grep $inte $DirPth/$F1/$F1.integrase.list.tmp | awk 'NR==1 {print $2}' )
            n1=$(grep -n "ORIGIN" $DirPth/$F1/$F1.re-annotated.gbk | awk -F':' '{print $1}')
            n2=$(( $n1 - 1))
            V1=$(echo '&&&&&misc_feature&&&&'$start'..'$end'==&&&&&&&&&&&&&&&&&&&&&/ugene_group='MGE'==&&&&&&&&&&&&&&&&&&&&&/ugene_name='MGE'==&&&&&&&&&&&&&&&&&&&&&/gene='integrases'==&&&&&&&&&&&&&&&&&&&&&/note='integrases'' | sed "s/ //g")
            sed -i -e ''$n2'a\'$V1'' $DirPth/$F1/$F1.re-annotated.gbk
        done
    fi

    ##-----------------------------------------------------------------------------
    ## plasmid replicons
    echo "running detectoin of plasmid replicons for $F1"
    cat results/25_plasmid_replicons"$s"/raw_files/$F1.plasmid.blast.tmp | awk -F'_' '{print $1}' | sort -u > $DirPth/$F1/$F1.p_rep.list
    for p_rep2 in $(cat $DirPth/$F1/$F1.p_rep.list); do
        start=$(grep $p_rep2 results/25_plasmid_replicons/raw_files/$F1.plasmid.blast.tmp | awk '{print $5}' | head -1 )
        end=$(grep $p_rep2 results/25_plasmid_replicons/raw_files/$F1.plasmid.blast.tmp | awk '{print $6}' | head -1 )
        n1=$(grep -n "ORIGIN" $DirPth/$F1/$F1.re-annotated.gbk | awk -F':' '{print $1}')
        n2=$(( $n1 - 1))
        V1=$(echo '&&&&&misc_feature&&&&'$start'..'$end'==&&&&&&&&&&&&&&&&&&&&&/ugene_group='misc_feature'==&&&&&&&&&&&&&&&&&&&&&/ugene_name='misc_feature'==&&&&&&&&&&&&&&&&&&&&&/gene='$p_rep2'==&&&&&&&&&&&&&&&&&&&&&/note='$p_rep2'' | sed "s/ //g")
        sed -i -e ''$n2'a\'$V1'' $DirPth/$F1/$F1.re-annotated.gbk
    done

        ##-----------------------------------------------------------------------------
    ## Virulence genes
    echo "running VG for $F1"
    file1=$(find results/21_VG_diamond_results*/tmp2/ -name $F1.VFDB-blast.csv | head -1)
    n=$(wc -l $file1 | awk '{print $1}')
    if [ $n -gt 1 ] ; then
        awk 'NR>1 {print $2}' $file1 > $DirPth/$F1/$F1.VFDB-blast.list
        for VG in $(cat $DirPth/$F1/$F1.VFDB-blast.list); do
            start=$(grep $VG results/08_annotation/raw_files/$F1/$F1.gbk.csv | awk -F',' 'NR==1 {print $3}' | sed 's/\"//g')
            end=$(grep $VG results/08_annotation/raw_files/$F1/$F1.gbk.csv | awk -F',' 'NR==1 {print $4}' | sed 's/\"//g')
            VG_gene=$(grep $VG $file1 | awk -F'\t' '{print $13}' )
            VG_note=$(grep $VG $file1 | awk -F'\t' '{print $14,$15,$16}' | sed "s/ /_/g")
            n1=$(grep -n "ORIGIN" $DirPth/$F1/$F1.re-annotated.gbk | awk -F':' '{print $1}')
            n2=$(( $n1 - 1))
            V1=$(echo '&&&&&misc_feature&&&&'$start'..'$end'==&&&&&&&&&&&&&&&&&&&&&/ugene_group='VG'==&&&&&&&&&&&&&&&&&&&&&/ugene_name='VG'==&&&&&&&&&&&&&&&&&&&&&/gene='$VG_gene'==&&&&&&&&&&&&&&&&&&&&&/note='$VG_note'' | sed "s/ //g")
            sed -i -e ''$n2'a\'$V1'' $DirPth/$F1/$F1.re-annotated.gbk
        done
    fi
    (rm $DirPth/$F1/$F1.VFDB-blast.list) > /dev/null 2>&1
    ##-----------------------------------------------------------------------------
    ## Abr genes
    echo "running Abr genes for $F1"
    file2=$(find results/22_CARD-AR_results*/results/ -name $F1.txt | head -1)
    n=$(wc -l $file2 | awk 'NR==1 {print $1}')
    if [ $n -gt 1 ] ; then
        awk -F' ' 'NR>1 {print $1}' $file2 > $DirPth/$F1/$F1.CARD-AR.list
        for Abr in $(cat $DirPth/$F1/$F1.CARD-AR.list); do
            start=$(grep $Abr results/08_annotation/raw_files/$F1/$F1.gbk.csv | awk -F',' 'NR==1 {print $3}' | sed 's/\"//g')
            end=$(grep $Abr results/08_annotation/raw_files/$F1/$F1.gbk.csv | awk -F',' 'NR==1 {print $4}' | sed 's/\"//g')
            Abr_gene=$(grep $Abr $file2 | awk -F'\t' 'NR==1 {print $9}' | sed "s/ /_/g")
            Abr_note=$(grep $Abr $file2 | awk -F'\t' 'NR==1 {print $15,$16}' | sed "s/ /_/g")
            n1=$(grep -n "ORIGIN" $DirPth/$F1/$F1.re-annotated.gbk | awk -F':' '{print $1}')
            n2=$(( $n1 - 1))
            V1=$(echo '&&&&&misc_feature&&&&'$start'..'$end'==&&&&&&&&&&&&&&&&&&&&&/ugene_group='Abr'==&&&&&&&&&&&&&&&&&&&&&/ugene_name='Abr'==&&&&&&&&&&&&&&&&&&&&&/gene='$Abr_gene'==&&&&&&&&&&&&&&&&&&&&&/note='$Abr_note'' | sed "s/ //g")
            sed -i -e ''$n2'a\'$V1'' $DirPth/$F1/$F1.re-annotated.gbk
        done
    fi
    (rm $DirPth/$F1/$F1.CARD-AR.list) > /dev/null 2>&1

##-----------------------------------------------------------------------------
    ## Genomic island

    echo "running GI for $F1"
    n=$(wc -l results/43_island_finder/raw_result_files/$F1.islands.gff | awk '{print $1}')
    if [ $n -gt 1 ] ; then
        awk 'NR>1 {print "Gen_Island"NR-1, $0}' results/43_island_finder/raw_result_files/$F1.islands.gff > $DirPth/$F1/$F1.islands.prefixed.gff 
        awk '{print $1}' $DirPth/$F1/$F1.islands.prefixed.gff > $DirPth/$F1/$F1.islands.prefixed.gff.list
        for GI in $(cat $DirPth/$F1/$F1.islands.prefixed.gff.list); do
            start=$(grep $GI $DirPth/$F1/$F1.islands.prefixed.gff | awk '{print $5}' )
            end=$(grep $GI $DirPth/$F1/$F1.islands.prefixed.gff | awk '{print $6}' )
            n1=$(grep -n "ORIGIN" $DirPth/$F1/$F1.re-annotated.gbk | awk -F':' '{print $1}')
            n2=$(( $n1 - 1))
            V1=$(echo '&&&&&misc_feature&&&&'$start'..'$end'==&&&&&&&&&&&&&&&&&&&&&/note="GI"==&&&&&&&&&&&&&&&&&&&&&/ugene_name="GI"==&&&&&&&&&&&&&&&&&&&&&/ugene_group="GI"' | sed "s/ //g")
            sed -i -e ''$n2'a\'$V1'' $DirPth/$F1/$F1.re-annotated.gbk
        done
    fi
    (rm $DirPth/$F1/$F1.islands.prefixed.gff ) > /dev/null 2>&1
    (rm $DirPth/$F1/$F1.islands.prefixed.gff.list) > /dev/null 2>&1



    ##-----------------------------------------------------------------------------
    sed -i 's/&&&&&&&&&&&&&&&&&&&&&/                     /g' $DirPth/$F1/$F1.re-annotated.gbk
    sed -i 's/==/\n/g' $DirPth/$F1/$F1.re-annotated.gbk
    sed -i 's/&&&&&misc_feature&&&&/     misc_feature    /g' $DirPth/$F1/$F1.re-annotated.gbk
    ##-----------------------------------------------------------------------------

done
###############################################################################

echo "finished.. step-08 re-annotation ---------------------------------------"

###############################################################################

exit

##-----------------------------------------------------------------------------
    ## plasmid contigs
    
    echo "running plasmid for $F1"
    grep "$F1" results/25_plasmid_finder_complete_MASH"$s"/results.csv | awk '{ s = ""; for (i = 13; i <= NF; i++) s = s $i " "; print s }' | tr " " "\n" | sed -r '/^\s*$/d' | awk -F'_' '$4>1000' | awk -F'_' '$6>5' > $DirPth/$F1/$F1.plasmid_contigs.list

    n=$(wc -l $DirPth/$F1/$F1.plasmid_contigs.list | awk '{print $1}')
    if [ $n -gt 1 ] ; then
        perl /home/swapnil/pipeline/tools/fastagrep.pl -f $DirPth/$F1/$F1.plasmid_contigs.list results/04_assembly/raw_files/$F1/$F1.fasta | awk '/^>/ {F = $1".fa"} {print > F}' 
        rename 's/\>//' *.fa
        for pc in $(cat $DirPth/$F1/$F1.plasmid_contigs.list);do
            start_end=$(blastn -query $pc.fa -subject results/08_annotation/raw_files/$F1/$F1.fna -max_target_seqs 1 -max_hsps 1 -evalue 1e-100 -outfmt "6 sstart send")
            start=$(echo $start_end | awk '{print $1}')
            end=$(echo $start_end | awk '{print $2}')
            n1=$(grep -n "ORIGIN" $DirPth/$F1/$F1.re-annotated.gbk | awk -F':' '{print $1}')
            n2=$(( $n1 - 1))
            V1=$(echo '&&&&&misc_feature&&&&'$start'..'$end'==&&&&&&&&&&&&&&&&&&&&&/ugene_group='pContig'==&&&&&&&&&&&&&&&&&&&&&/ugene_name='pContig'' | sed "s/ //g")
            sed -i -e ''$n2'a\'$V1'' $DirPth/$F1/$F1.re-annotated.gbk
        done
        rm *.fa
        #(rm results/08_annotation/00_re-annotated/$F1/$F1.plasmid_contigs.list) > /dev/null 2>&1
    fi

    