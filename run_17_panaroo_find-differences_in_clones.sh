#!/bin/bash
###############################################################################
#9 WGS-alignment-mugsy
# add -z option to keep the intermediate files
###############################################################################
## note: 
    #echo "provide the list of clone_id and respective isolates_id\n the list should be serial starting from 1 "
    #echo "the tab seperated file (clones.txt) should contain"
    #echo "-------------------------------------------------------------------------"
    #echo "1 strain1"
    #echo "1 strain2"
    #echo "2 strain5"
    #echo "2 strain9"
    #echo "3 strain10"
    #echo "4 strain15"
    #echo "4 strain19"
    #echo "-------------------------------------------------------------------------"
###############################################################################
## file and directory preparation step 

    (mkdir results/17_roary_differences_in_clones) > /dev/null 2>&1

    path1=results/17_roary_differences_in_clones
###############################################################################
## seperate genes
    #echo "seperating genes"
    cat clones.txt | awk '{print $1}' | sort -nu > $path1/clones.list
    total_No_clones=$(cat clones.txt | awk '{print $1}' | sort -u | wc -l)

    clone_no=1
    while : ; do
        if [ -f $path1/$clone_no/"Clone_"$clone_no.csv ]; then
            echo "gene seperation already finished for $clone_no"
            clone_no=$[ $clone_no + 1 ]
            else
            # echo "gene seperation running for clone $clone_no"
            if [ "$clone_no" -gt "$total_No_clones" ]; then
                #echo "breaking"
                break
                else
                if [ ! -d $path1/$clone_no ] ; then
                    echo "running gene seperation for Clone $clone_no ========================================="
                    (mkdir $path1/$clone_no) > /dev/null 2>&1
                    (mkdir $path1/$clone_no/tmp) > /dev/null 2>&1
                    (mkdir $path1/$clone_no/gff) > /dev/null 2>&1 
                    cat clones.txt | awk '$1 == '$clone_no' {print $2}' | sed 's/ /\n/g' > clones.$clone_no.txt.tmp
                    total_No_clones2=$( wc -l clones.$clone_no.txt.tmp | awk '{print $1}' )
                    pan_genome=$(( $total_No_clones2 - 1 ))
                    #echo "$total_No_clones2 $pan_genome pan_genome"

                    for isolate in $(cat clones.$clone_no.txt.tmp) ; do
                        echo "copying $isolate"
                        cp results/08_annotation/raw_files/$isolate/$isolate.gff $path1/$clone_no/gff/
                        #cp /home/swapnil/p_Entb_Germany/results/17_roary_165/gff/$isolate.gff $path1/$clone_no/gff/
                    done

                    ## run roary
                    # (roary -p 8 -f $path1/$clone_no/roary_results -i 70 $path1/$clone_no/gff/*.gff) > /dev/null 2>&1
                    ## run panaroo
                    echo "running panaroo for clone $clone_no"
                    (mkdir $path1/$clone_no/roary_results) > /dev/null 2>&1
                    panaroo \
                    -i $path1/$clone_no/gff/*.gff \
                    -o $path1/$clone_no/roary_results \
                    -t 8 \
                    -c 0.70 \
                    --len_dif_percent 0.7 \
                    --quiet \
                    --merge_paralogs \
                    --clean-mode strict

                    cat $path1/$clone_no/roary_results/gene_presence_absence_roary.csv | awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", "", $i) } 1' | sed 's/"//g' | awk -F',' '$4<="'$pan_genome'" {print $0}' > $path1/$clone_no/tmp/gene_presence_absence.csv.tmp
                    head -1 $path1/$clone_no/roary_results/gene_presence_absence_roary.csv | sed 's/"//g' | sed 's/ //g' > $path1/$clone_no/tmp/head.tmp
                    cat $path1/$clone_no/tmp/head.tmp $path1/$clone_no/tmp/gene_presence_absence.csv.tmp > $path1/$clone_no/tmp/gene_presence_absence.csv.tmp2
                    csvcut --columns "Gene,Annotation" $path1/$clone_no/tmp/gene_presence_absence.csv.tmp2 > $path1/$clone_no/tmp/gene_annotation.tmp

                    for isolate in $(cat clones.$clone_no.txt.tmp) ; do
                        csvcut --columns "$isolate" $path1/$clone_no/tmp/gene_presence_absence.csv.tmp2 > $path1/$clone_no/tmp/$isolate.txt
                    done
                    
                    paste $path1/$clone_no/tmp/gene_annotation.tmp $path1/$clone_no/tmp/*.txt > $path1/$clone_no/"Clone_"$clone_no.csv

                    (rm clones.$clone_no.txt.tmp) > /dev/null 2>&1
                    echo "Gene seperation finished for $clone_no =============================================="
                    clone_no=$[ $clone_no + 1 ]
                    else
                    echo "gene seperation step for clone $clone_no already finished"
                    clone_no=$[ $clone_no + 1 ]
                fi
            fi
        fi
    done
###############################################################################
## Plasmid
    for clone_no in $(cat $path1/clones.list); do
        (rm $path1/$clone_no/tmp/gene_name.panaroo_refound.tmp) > /dev/null 2>&1
        if [ -f $path1/$clone_no/Clone_$clone_no.plasmid.csv ] ; then 
            echo "plasmid analysis for clone $clone_no already finished"
            else 
            echo "started for clone $clone_no -----------------------------------------------"
            awk -F',' 'NR> 1 {print $1}' $path1/$clone_no/"Clone_"$clone_no.csv | awk -F'\t' '{print $1}' > $path1/$clone_no/tmp/"Clone_"$clone_no.genes.csv

            ## in panaroo duplicated/truncated genes are labelled as refound and need to remove
            (rm $path1/$clone_no/tmp/"Clone_"$clone_no.genes.2.csv) > /dev/null 2>&1
            for gene in $(cat $path1/$clone_no/tmp/"Clone_"$clone_no.genes.csv); do
                gene_name=$(grep -w "$gene" $path1/$clone_no/"Clone_"$clone_no.csv | awk -F'\t' '{$1=""; print $0}' | sed 's/""//g' | awk '{print $1}' | awk -F';' '{print $1}')
                gene_name2=$( echo "$gene_name" | grep 'refound' )
                if [ -z $gene_name2 ] ; then
                    #echo gene_name $gene_name
                    echo $gene_name >> $path1/$clone_no/tmp/"Clone_"$clone_no.genes.2.csv
                    else
                    cat clones.txt | awk '$1 == '$clone_no' {print $2}' | sed 's/ /\n/g' > clones.$clone_no.txt.tmp
                    (rm $path1/$clone_no/tmp/gene_name.panaroo_refound.tmp) > /dev/null 2>&1
                    for isolate2 in $(cat clones.$clone_no.txt.tmp ); do
                        V1=$(grep -w "$gene" $path1/$clone_no/Clone_$clone_no.csv | nawk -F "$isolate2" '{print $2}' | awk '{print $1}' | awk -F';' '{print $1}')
                        if [ ! -z "$V1" ] ; then 
                            V2=$(echo "$isolate2""$V1")
                            #echo V2 $V2
                            echo "$V2" >> $path1/$clone_no/tmp/gene_name.panaroo_refound.tmp
                        fi
                    done
                    V3=$(head -1 $path1/$clone_no/tmp/gene_name.panaroo_refound.tmp)
                    echo $V3 >> $path1/$clone_no/tmp/"Clone_"$clone_no.genes.2.csv
                    rm $path1/$clone_no/tmp/gene_name.panaroo_refound.tmp
                fi
            done
            ## BLAST for plasmid
            echo "..... running plasmid-BLAST"
            (rm $path1/$clone_no/tmp/plasmid_genes.txt) > /dev/null 2>&1
            for gene_2 in $(cat $path1/$clone_no/tmp/"Clone_"$clone_no.genes.2.csv) ; do
                isolate=$(echo $gene_2 | sed 's/......$//' )
                faidx results/08_annotation/raw_files/$isolate/$isolate.ffn $gene_2 > $path1/$clone_no/tmp/$gene_2.fasta
                blast_results=$(blastn -db /media/swapnil/share/databases/plasmid_06022020_brooks_curated/all_plasmids.fna -query $path1/$clone_no/tmp/$gene_2.fasta -max_target_seqs 5 -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk '$14>85 && $15>70' | awk 'NR==1 {print $2}')
                    if [ -z "$blast_results" ] ; then
                    echo -e "$gene_2 \t 0" >> $path1/$clone_no/tmp/plasmid_genes.txt
                    else
                    echo -e "$gene_2 \t plasmid_gene" >> $path1/$clone_no/tmp/plasmid_genes.txt
                    fi
            done
            sed  -i '1i gene-BLASTed\tplasmid-gene' $path1/$clone_no/tmp/plasmid_genes.txt
            
            paste $path1/$clone_no/Clone_$clone_no.csv $path1/$clone_no/tmp/plasmid_genes.txt > $path1/$clone_no/Clone_$clone_no.plasmid.csv
        fi
    done
###############################################################################
## prophage
    for clone_no in $(cat $path1/clones.list); do
        if [ -f $path1/$clone_no/Clone_$clone_no.plasmid.prophage.csv ]; then
            echo "prophage-BLAST for clone $clone_no already finished"
            else
            ## BLAST for prophage
            echo "..... running prophage-BLAST for clone $clone_no"
            (rm $path1/$clone_no/tmp/prophage_genes.txt) > /dev/null 2>&1
            for gene_2 in $(cat $path1/$clone_no/tmp/"Clone_"$clone_no.genes.2.csv) ; do
            already_a_plasmid_gene=$(grep $gene_2 $path1/$clone_no/Clone_$clone_no.plasmid.csv | awk '{print $NF}')
                if [ "$already_a_plasmid_gene" == "plasmid_gene" ]; then
                    echo -e "$gene_2 \t 0" >> $path1/$clone_no/tmp/prophage_genes.txt
                    else
                    blast_results=$(blastn -db /media/swapnil/share/databases/prophage_06022020/all_prophages.fasta -query $path1/$clone_no/tmp/$gene_2.fasta -max_target_seqs 5 -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk '$14>85 && $15>70' | awk 'NR==1 {print $2}')
                        if [ -z "$blast_results" ] ; then
                            echo -e "$gene_2 \t 0" >> $path1/$clone_no/tmp/prophage_genes.txt
                            else
                            echo -e "$gene_2 \t prophage_gene" >> $path1/$clone_no/tmp/prophage_genes.txt
                        fi
                fi
            done
            sed  -i '1i gene-BLASTed\tprophage-gene' $path1/$clone_no/tmp/prophage_genes.txt
            paste $path1/$clone_no/Clone_$clone_no.plasmid.csv $path1/$clone_no/tmp/prophage_genes.txt > $path1/$clone_no/Clone_$clone_no.plasmid.prophage.csv
        fi
    done
###############################################################################
## find gap_effect
    ## re-format fasta (upstram-downsrtream region extraction step need reformated fasta)
    cat clones.txt | awk '{print $2}' | sort -u > reformat_fasta_list.txt
    ( mkdir $path1/tmp )> /dev/null 2>&1
    ( mkdir $path1/tmp/fasta )> /dev/null 2>&1

    for fasta in $(cat reformat_fasta_list.txt); do
        #echo "reformatting $fasta fasta"
        cp results/04_assembly/all_fasta/$fasta.fasta $path1/tmp/fasta/$fasta.fasta
        (python /home/swapnil/tools/read-cleaning-format-conversion/KSU_bioinfo_lab/fasta-o-matic/fasta_o_matic.py -o $path1/tmp/fasta/ -f $path1/tmp/fasta/$fasta.fasta)> /dev/null 2>&1
        (mv $path1/tmp/fasta/"$fasta"_wrap.fasta $path1/tmp/fasta/$fasta.fasta)> /dev/null 2>&1
    done
    ( rm reformat_fasta_list.txt )> /dev/null 2>&1
    #echo "reformatting fasta step finished"

    for clone_no in $(cat $path1/clones.list); do
        if [ -f $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.csv ] ; then
            echo "gap effect for $clone_no already finished"
            else
            echo ".... running gap effect for $clone_no"
            ( rm $path1/$clone_no/tmp/gap_effect_genes.txt )> /dev/null 2>&1
            ( mkdir $path1/$clone_no/tmp )> /dev/null 2>&1
            ( mkdir $path1/$clone_no/fasta_extracted )> /dev/null 2>&1
            ( mkdir $path1/$clone_no/fasta_extracted/discarded )> /dev/null 2>&1
            
            for gene_2 in $(cat $path1/$clone_no/tmp/"Clone_"$clone_no.genes.2.csv) ; do
                F1=$( echo $gene_2 | awk -F'_' '{print $1}')
                #echo "blasting $F1 $gene_2"
                blastn -subject $path1/tmp/fasta/$F1.fasta -query $path1/$clone_no/tmp/$gene_2.fasta -out $path1/tmp/fasta/$F1.$gene_2.custom-gene-blast.tmp -max_target_seqs 1 -max_hsps 1 -evalue 1e-10 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
                #if [ -f $path1/tmp/fasta/$F1.$gene_2.custom-gene-blast.tmp ] ; then echo "blastn succesful" ; else echo "blastn failed" ; fi ;
                if [ -s "$path1/tmp/fasta/$F1.$gene_2.custom-gene-blast.tmp" ]; then
                    start=$(awk 'NR==1 {print $3}' $path1/tmp/fasta/$F1.$gene_2.custom-gene-blast.tmp)
                    end=$(awk 'NR==1 {print $4}' $path1/tmp/fasta/$F1.$gene_2.custom-gene-blast.tmp)
                    if [ "$start" -gt "$end" ]; then
                        #echo "$F1 $gene_2 reverse"
                        Forw=$(( $end - 100 ))
                        Revr=$(( $start + 100 ))
                        (rm $path1/tmp/fasta/$F1.fasta.fai)> /dev/null 2>&1
                        echo -e "$F1\t$Forw\t$Revr\t"$F1"_region" > $path1/tmp/fasta/$F1.$gene_2.bed
                        (bedtools getfasta -fi $path1/tmp/fasta/$F1.fasta -bed $path1/tmp/fasta/$F1.$gene_2.bed -name > $path1/$clone_no/fasta_extracted/"$gene_2"_region.RC.fasta )> /dev/null 2>&1
                        sed '1d' $path1/$clone_no/fasta_extracted/"$gene_2"_region.RC.fasta | perl -pe 'chomp;tr/ACGTNacgtn/TGCANtgcan/;$_=reverse."\n"' | sed '1i >'$F1'' > $path1/$clone_no/fasta_extracted/"$gene_2"_region.fasta
                        mv $path1/$clone_no/fasta_extracted/"$gene_2"_region.RC.fasta $path1/$clone_no/fasta_extracted/discarded/
                        else
                        #echo "$F1 $gene_2 normal"
                        Forw=$(( $start - 100 ))
                        Revr=$(( $end + 100 ))
                        (rm $path1/tmp/fasta/$F1.fasta.fai)> /dev/null 2>&1
                        echo -e "$F1\t$Forw\t$Revr\t"$F1"_region" > $path1/tmp/fasta/$F1.$gene_2.bed
                        (bedtools getfasta -fi $path1/tmp/fasta/$F1.fasta -bed $path1/tmp/fasta/$F1.$gene_2.bed -name > $path1/$clone_no/fasta_extracted/"$gene_2"_region.fasta )> /dev/null 2>&1
                    fi
                    (rm $path1/tmp/$F1.custom-gene-blast.tmp)> /dev/null 2>&1
                fi
                #echo "gene upstram-downstream extracted"

                gap=$(cat $path1/$clone_no/fasta_extracted/"$gene_2"_region.fasta | grep -c NNNNNNNNNN )
                if [ "$gap" -gt "0" ] ; then
                    echo -e "$gene_2 \t gap_effect" >> $path1/$clone_no/tmp/gap_effect_genes.txt
                    else
                    echo -e "$gene_2 \t 0" >> $path1/$clone_no/tmp/gap_effect_genes.txt
                fi
            done
        sed  -i '1i gene-BLASTed\tgap_effect' $path1/$clone_no/tmp/gap_effect_genes.txt
        paste $path1/$clone_no/Clone_$clone_no.plasmid.prophage.csv $path1/$clone_no/tmp/gap_effect_genes.txt > $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.csv

        ( rm -rf $path1/$clone_no/fasta_extracted )>/dev/null 2>&1
        ( rm -rf $path1/$clone_no/fasta_extracted/discarded )>/dev/null 2>&1
        fi
    done
    ( rm -rf $path1/tmp ) >/dev/null 2>&1
    ( rm -rf $path1/tmp/fasta )>/dev/null 2>&1
    ( rm -rf $path1/tmp/fasta_extracted ) >/dev/null 2>&1
###############################################################################
## Key genes
    for clone_no in $(cat $path1/clones.list); do
        if [ -f $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.csv ]; then
            echo "keygenes for clone $clone_no already finished"
            else 
            echo "..... running keygenes for clone $clone_no"
            (rm $path1/$clone_no/tmp/*.keygenes.txt) > /dev/null 2>&1
            (rm $path1/$clone_no/tmp/all.keygenes.txt) > /dev/null 2>&1
            ##
            grep "tRNA" $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.csv | awk -F',' '{print $0, "tRNA"}' > $path1/$clone_no/tmp/tRNA.keygenes.txt
            grep "transposase" $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.csv | awk -F',' '{print $0, "mobile-genetic-element"}' > $path1/$clone_no/tmp/transposase.keygenes.txt
            ##
            grep "Prophage integrase" $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.csv | awk -F',' '{print $0, "prophage_gene"}' > $path1/$clone_no/tmp/prophage_integrase.keygenes.txt
            grep "recombinase" $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.csv | awk -F',' '{print $0, "prophage_gene"}' > $path1/$clone_no/tmp/prophage_recombinase.keygenes.txt
            ## plasmid
            grep "replication protein" $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.csv | awk -F',' '{print $0, "plasmid_gene"}' > $path1/$clone_no/tmp/replication_protein.keygenes.txt
            grep "plasmid" $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.csv | awk -F',' '{print $0, "plasmid_gene"}' > $path1/$clone_no/tmp/plasmid.keygenes.txt

            cat $path1/$clone_no/tmp/*.keygenes.txt | tr ' ' '\t' > $path1/$clone_no/tmp/all.keygenes.txt

            (rm $path1/$clone_no/tmp/keygenes.txt) > /dev/null 2>&1
            for gene_2 in $(cat $path1/$clone_no/tmp/"Clone_"$clone_no.genes.2.csv) ; do
                keygene=$(grep $gene_2 $path1/$clone_no/tmp/all.keygenes.txt | awk -F'\t' '{print $NF}' | tr '\r\n' '&' )
                keygene=$(grep $gene_2 $path1/$clone_no/tmp/all.keygenes.txt | head -1 | awk -F'\t' '{print $NF}' )
                if [ -z "$keygene" ] ; then ## if empty
                    echo -e "$gene_2 \t 0" >> $path1/$clone_no/tmp/keygenes.txt
                    else
                    echo -e "$gene_2 \t $keygene" >> $path1/$clone_no/tmp/keygenes.txt
                fi
            done
            sed -i '1i gene-BLASTed\tKeygene' $path1/$clone_no/tmp/keygenes.txt
            paste $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.csv $path1/$clone_no/tmp/keygenes.txt > $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.csv
        fi
    done
###############################################################################
## summarise all genes
    (rm $path1/results.txt) > /dev/null 2>&1
    (rm $path1/CloneIsolate.results.txt) > /dev/null 2>&1
    for clone_no in $(cat $path1/clones.list); do
        awk '{print $(NF-6), $(NF-4), $(NF-2), $(NF)}' $path1/$clone_no/Clone_$clone_no.plasmid.prophage.gap_effect.keygenes.csv | sed 's/0 //g' | sed 's/^0/,/g' | sed 's/0$//g' | sed 's/,//g' | sed 's/ / \& /g' | sed 's/ \& $//g' | sed '1 s/^.*$/transgenetic-element/' > $path1/$clone_no/tmp/allgenes.csv.tmp
        (rm $path1/$clone_no/tmp/allgenes.csv) > /dev/null 2>&1
        while read duplicates; do 
        echo "$duplicates" | awk '{for (i=1;i<=NF;i++) if (!a[$i]++) printf("%s%s",$i,FS)}{printf("\n")}' | sed 's/ & $//g' >> $path1/$clone_no/tmp/allgenes.csv
        done < $path1/$clone_no/tmp/allgenes.csv.tmp

        paste $path1/$clone_no/Clone_$clone_no.csv $path1/$clone_no/tmp/allgenes.csv > $path1/$clone_no/tmp/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.csv

        sed '/gap_effect/d' $path1/$clone_no/tmp/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.csv > $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv
        
        core_gene=$(cat $path1/$clone_no/roary_results/summary_statistics.txt | head -1 | awk '{print $NF}')
        total_genes=$(awk 'NR>1 {print $NF}' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv | wc -l)
        pan_genes=$(( $core_gene + $total_genes ))
        transgenetic_genes=$(awk -F'\t' 'NR>1 {print $NF}' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv | sed '/^\s*$/d' | wc -l)
        true_genes="$(( $total_genes - $transgenetic_genes ))"        #percentage_of_transgentic_genes=$(echo 'scale=1;'$transgenetic_genes'*'100'/'$total_genes'' | bc)

        plasmid_gene=$(grep -c 'plasmid_gene' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv)
        prophage_gene=$(grep -c 'prophage_gene' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv)
        mobile_genetic_element=$(grep -c 'mobile-genetic-element' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv)

        echo "Clone_$clone_no, core/pan_genes: $core_gene/$pan_genes, genes difference = $total_genes (plasmid_genes = $plasmid_gene + prophage_gene = $prophage_gene + mobile-genetic-element = $mobile_genetic_element), transgenetic genes = $transgenetic_genes, true-gene-difference-upto = $true_genes" > $path1/$clone_no/result.txt
        echo "Clone_$clone_no, core/pan_genes: $core_gene/$pan_genes, genes difference = $total_genes (plasmid_genes = $plasmid_gene + prophage_gene = $prophage_gene + mobile-genetic-element = $mobile_genetic_element), transgenetic genes = $transgenetic_genes, true-gene-difference-upto = $true_genes" >> $path1/results.txt


        ## results for each isolate of each clone
        cat clones.txt | awk '{print $1}' | sort -nu > clones.list

        awk '$1 == '$clone_no' {print $2}' clones.txt > $path1/$clone_no/tmp/list.isolate.txt
        grep 'prophage_gene' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv >  $path1/$clone_no/tmp/final.prophage_gene.csv
        grep 'plasmid_gene' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv >  $path1/$clone_no/tmp/final.plasmid_gene.csv
        grep 'mobile-genetic-element' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv >  $path1/$clone_no/tmp/final.mge_gene.csv
        grep 'tRNA' $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv >  $path1/$clone_no/tmp/final.tRNA_gene.csv


        for isolate in $(cat $path1/$clone_no/tmp/list.isolate.txt); do
            #isolate_pan_gene=$(cat $path1/$clone_no/tmp/gene_presence_absence.csv.tmp | grep $isolate | wc -l )
            gap_effect_genes=$(cat $path1/$clone_no/tmp/gap_effect_genes.txt | grep $isolate | grep gap_effect | wc -l )
            isolate_plasmid=$(cat $path1/$clone_no/tmp/final.plasmid_gene.csv | grep $isolate | wc -l )
            isolate_prophage=$(cat $path1/$clone_no/tmp/final.prophage_gene.csv | grep $isolate | wc -l)    
            isolate_mge=$(cat $path1/$clone_no/tmp/final.mge_gene.csv | grep $isolate | wc -l)
            isolate_tRNA=$(cat $path1/$clone_no/tmp/final.tRNA_gene.csv | grep $isolate | wc -l)
            isolate_mge2=$(( $isolate_mge + $isolate_tRNA ))
            other_genes=$(cat $path1/$clone_no/Clone_"$clone_no".plasmid.prophage.gap_effect.keygenes.allgenes.gap_effect_removed.csv | grep $isolate | awk -F'\t' '{print $NF}' | grep -c "^$")
            isolate_pan_gene=$(($gap_effect_genes + $isolate_plasmid + $isolate_prophage + $isolate_mge2 + $other_genes))

            echo $clone_no $isolate $core_gene $isolate_pan_gene $gap_effect_genes $isolate_plasmid $isolate_prophage $isolate_mge2 $other_genes > $path1/$clone_no/CloneIsolate.result.txt
            echo $clone_no $isolate $core_gene $isolate_pan_gene $gap_effect_genes $isolate_plasmid $isolate_prophage $isolate_mge2 $other_genes >> $path1/CloneIsolate.results.txt
        done
    done

    sed -i '1i clone_no isolate core_genes pan_genes truncated_due_to_assembly plasmid prophage MGE other_genes' $path1/CloneIsolate.results.txt

###############################################################################
