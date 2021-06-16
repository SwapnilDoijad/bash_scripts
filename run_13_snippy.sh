#!/bin/bash
###############################################################################
#13 snippy
## FOR GENOMEWIDE: script is only for core SNVs. filter 0 SNVs (0 is producing wrong average). need to work on averages 
###############################################################################
## preliminary file preparation and directory creation
    echo "started.... step-13 snippy ----------------------------------------------"

    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "want to suffix files? type "y" else press enter to continue"
    read answer
    if [ "$answer" = "y" ]; then
    s=$(echo "_$l")
    fi

    echo "choose reference, provide full path"
    ls results/00_ref/*.gbk | sed 's/ /\n/g'
    read ref

    echo "If wish to use (filtered-)reads type 'y' else press ENTER to continue"
    read answer_2
    if [ "$answer_2" == "y" ] ; then
    echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany) provide full path"
    read path
    fi

    #echo "Wish to run the detailed SNP analysis? If yes type y or press Enter"
    #read answer2    

    ( mkdir results/13_snippy"$s" ) > /dev/null 2>&1
    ( mkdir results/13_snippy"$s"/tmp ) > /dev/null 2>&1
###############################################################################
## run snippy program

    for F1 in $(cat $list); do
        if [ -d results/13_snippy"$s"/$F1 ] ; then
        echo "$F1 is already done" 
        else
            echo "running snippy for $F1"
            /home/swapnil/tools/snippy/bin/snippy --cpus 16 --quiet --prefix $F1 --force --outdir  results/13_snippy"$s"/$F1 --ref $ref --pe1 $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz --pe2 $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz #) > /dev/null 2>&1
            #fasta
            #(snippy --cpus 16 --quiet --prefix $F1 --force --outdir  results/13_snippy"$s"/$F1 --ref $ref --ctgs results/04_assembly/all_fasta/$F1.fasta) > /dev/null 2>&1
        fi
    done

    (rm list.tmp ) > /dev/null 2>&1
    for F1 in $(cat $list); do
        echo results/13_snippy"$s"/$F1 >> list.tmp
    done

    /home/swapnil/tools/snippy/bin/snippy-clean_full_aln results/13_snippy"$s"/core.full.aln > results/13_snippy"$s"/clean.full.aln
    #run_gubbins.py -p gubbins results/13_snippy"$s"/clean.full.aln
    snp-sites -b -c -o results/13_snippy"$s"/clean.core.aln results/13_snippy"$s"/clean.full.aln
    #FastTree -gtr -nt results/13_snippy"$s"/clean.core.aln > results/13_snippy"$s"clean.core.aln.tree

    /home/swapnil/tools/snippy/bin/snippy-core --prefix core --ref $ref $(cat list.tmp)
    mv core.* results/13_snippy"$s"/
    (rm list.tmp)> /dev/null 2>&1

    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate myenv
    snp-dists -q results/13_snippy"$s"/core.full.aln > results/13_snippy"$s"/clean.full.aln.matrix.tsv
    conda deactivate

exit
###############################################################################
## detailed analysis
    if [ "$answer2" == "y" ] ; then

        ##extract ref LOCUS_TAG
        ref=$(ls results/00_ref/ref.*.gbk )
        (rm "$ref".csv) > /dev/null 2>&1
        (ugene --task=/home/swapnil/pipeline/tools/gbk2csv-Ugene.uwl --in=$ref --out="$ref".csv --format=csv) > /dev/null 2>&1 ;

    
        for F1 in $(cat $list); do
            grep "CDS" results/13_snippy"$s"/And1463/And1463.tab | awk -F'\t' 'FNR>1 {print $12}' | sort -u -k1,1 | sed '/^\s*$/d' > results/13_snippy"$s"/tmp/$F1.locus.tmp
        done

        for F1 in $(cat $list);do
            (mkdir results/13_snippy"$s"/$F1/tmp)> /dev/null 2>&1
            (rm results/13_snippy"$s"/tmp/01_SNPs_per_genes.tmp)> /dev/null 2>&1
            echo "counting SNPs in the genes of $F1"
            echo "gene length SNPs SNPs% stop_gained dN/dS product" > results/13_snippy"$s"/$F1/$F1.13_snippy.SNPs-statistics.txt
            for F2 in $(cat results/13_snippy"$s"/tmp/$F1.locus.tmp); do
                grep "$F2" results/13_snippy"$s"/$F1/$F1.tab > results/13_snippy"$s"/$F1/tmp/$F2.tab
                A0=$(grep "$F2" "$ref".csv | awk -F',' 'FNR==1 {print $5}' | sed "s/\"//g")
                A1=$(grep -c "$F2" results/13_snippy"$s"/$F1/tmp/$F2.tab) #SNPs per gene
                A2=$(awk "BEGIN {print 100 * $A1 / $A0 }" | awk '{printf "%.2f\n", $1, $2}'); if [ -z "$A2" ] ; then A2=$(echo "0") ; fi > /dev/null 2>&1
                A3=$(grep "$F2" "$ref".csv | awk -F',' 'FNR==1 {print $9}' | sed "s/\"//g" | sed "s/ /_/g")

                G0=$(grep -c 'initiator_codon_variant' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G0" ] ; then G0=$(echo "0") ; fi
                G1=$(grep -c 'initiator_codon_variant&non_canonical_start_codon' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G1" ] ; then G1=$(echo "0") ; fi
                G2=$(grep -c 'missense_variant' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G2" ] ; then G2=$(echo "0") ; fi
                G3=$(grep -c 'non_coding_transcript_variant' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G3" ] ; then G3=$(echo "0") ; fi
                G4=$(grep -c 'splice_region_variant&stop_retained_variant' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G4" ] ; then G4=$(echo "0") ; fi
                G5=$(grep -c 'start_lost' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G5" ] ; then G5=$(echo "0") ; fi
                G6=$(grep -c 'stop_gained' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G6" ] ; then G6=$(echo "0") ; fi
                G7=$(grep -c 'stop_lost&splice_region_variant' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G7" ] ; then G7=$(echo "0") ; fi
                G8=$(grep -c 'synonymous_variant' results/13_snippy"$s"/$F1/tmp/$F2.tab) ; if [ -z "$G8" ] ; then G8=$(echo "1") ; fi		#Synonymous
                G9=$(($G0+$G1+$G2+$G3+$G4+$G5+$G6+$G7))					#Sum of all non-synonymous SNPs
                G10=$(awk "BEGIN {print $G9 / $G8 }" | awk '{printf "%.2f\n", $1, $2}') ; if [ -z "$G10" ] ; then G10=$(echo "0") ; fi > /dev/null 2>&1	#dN/dS non-synonymous to synonnymous ratio
                echo $F2 $A0 $A1 $A2 $G6 $G10 $A3 >> results/13_snippy"$s"/$F1/$F1.13_snippy.SNPs-statistics.txt
            done
        done

        ###############################################################################
        ## genomewide

        (rm results/13_snippy"$s"/tmp/ref.gene.list)> /dev/null 2>&1
        cat $ref.csv | head -1 > $ref.CDS.csv
        cat $ref.csv | grep "CDS" >> $ref.CDS.csv
        cat $ref.CDS.csv | tr -s ' ' ',' | csvcut -c locus_tag | sed "s/\"//g" | sed -r '/^\s*$/d' >> results/13_snippy"$s"/tmp/ref.gene.list
        sed -i '/locus_tag/d' results/13_snippy"$s"/tmp/ref.gene.list # probelm if I combine with above line

        awk 'FNR>6 {print $2}' results/13_snippy"$s"/core.vcf > results/13_snippy"$s"/tmp/total_SNVs.list
        (rm results/13_snippy"$s"/tmp/total_SNVs.list.gene_found)> /dev/null 2>&1
        cp results/13_snippy"$s"/tmp/total_SNVs.list results/13_snippy"$s"/tmp/total_SNVs.list.tmp

        for gene in $(cat results/13_snippy"$s"/tmp/ref.gene.list); do
            start=$(grep "$gene" "$ref".CDS.csv | awk -F',' 'FNR==1 {print $3}' | sed "s/\"//g")
            end=$(grep "$gene" "$ref".CDS.csv | awk -F',' 'FNR==1 {print $4}' | sed "s/\"//g")

                for SNV in $(cat results/13_snippy"$s"/tmp/total_SNVs.list.tmp); do
                    if [[ "$SNV" -ge "$start" && "$SNV" -le "$end" ]] ; then
                    V1=$(grep "$gene" "$ref".CDS.csv | awk -F',' '{print $8}' | sed "s/\"//g")
                    echo "$SNV $V1"
                    echo "start"$SNV"end" $V1 >> results/13_snippy"$s"/tmp/total_SNVs.list.gene_found
                    fi
                done
        done

        echo "running post-processing"
        awk '{print $1}' results/13_snippy"$s"/tmp/total_SNVs.list.gene_found | sed "s/start//g" | sed "s/end//g" > results/13_snippy"$s"/tmp/total_SNVs.list.gene_found.2
        cat results/13_snippy"$s"/tmp/total_SNVs.list.gene_found.2 results/13_snippy"$s"/tmp/total_SNVs.list.tmp | awk '{!seen[$0]++};END{for(i in seen) if(seen[i]==1)print i}' | sort -n | awk 'NF{print "start"$1"end" " intergenic"}' > results/13_snippy"$s"/tmp/total_SNVs.list.not_found
        cat results/13_snippy"$s"/tmp/total_SNVs.list.gene_found results/13_snippy"$s"/tmp/total_SNVs.list.not_found > results/13_snippy"$s"/tmp/total_SNVs.list.all

        echo "running post-processing: creating total_SNVs.list.all.final"
        (rm results/13_snippy"$s"/tmp/total_SNVs.list.all.final)> /dev/null 2>&1
        awk 'NF{print "start"$0"end"}' results/13_snippy"$s"/tmp/total_SNVs.list > results/13_snippy"$s"/tmp/total_SNVs.list.2
        for F1 in $(cat results/13_snippy"$s"/tmp/total_SNVs.list.2); do
            grep "$F1" results/13_snippy"$s"/tmp/total_SNVs.list.all | awk 'FNR==1 {print $0}' >> results/13_snippy"$s"/tmp/total_SNVs.list.all.final
        done
        sed -i "s/start//g" results/13_snippy"$s"/tmp/total_SNVs.list.all.final
        sed -i "s/end//g" results/13_snippy"$s"/tmp/total_SNVs.list.all.final

        awk 'FNR>6 {print $0}' results/13_snippy"$s"/core.vcf > results/13_snippy"$s"/tmp/core.vcf.tmp

        paste results/13_snippy"$s"/tmp/total_SNVs.list.all.final results/13_snippy"$s"/tmp/core.vcf.tmp > results/13_snippy"$s"/core.vcf.annotated

        ##---------------------------------------------------
        cat $ref.CDS.csv | tr -s ' ' ',' | csvcut -c product | sed "s/\"//g" | sed -r '/^\s*$/d' | sed "s/,/_/g" > results/13_snippy"$s"/tmp/ref.gene_product.list ;
        sed -i '1d' results/13_snippy"$s"/tmp/ref.gene_product.list ;
        paste results/13_snippy"$s"/tmp/ref.gene.list results/13_snippy"$s"/tmp/ref.gene_product.list > results/13_snippy"$s"/tmp/ref.gene_product.list.2
        total_number_of_column=$(awk 'FNR>6 {print NF; exit}' results/13_snippy"$s"/core.vcf)
        total_number_of_strains=$(( $total_number_of_column - 9 ))
        awk '{for(i=12;i<=NF;i++) t+=$i; print $2, $1, t; t=0}' results/13_snippy"$s"/core.vcf.annotated > results/13_snippy"$s"/tmp/core.vcf.2

        (rm results/13_snippy"$s"/tmp/gene_locus.Poly_sites.tab) > /dev/null 2>&1
        for gene_locus in $(cat results/13_snippy"$s"/tmp/ref.gene.list); do
            length=$(grep "$gene_locus" "$ref".CDS.csv | awk -F',' '{print $5}' | sed "s/\"//g") ;
            gene_product=$(grep "$gene_locus" results/13_snippy"$s"/tmp/ref.gene_product.list.2 | awk '{print $2}' | sed "s/\"//g");
            Number_of_Polymorphic_sites=$(grep -c "$gene_locus" results/13_snippy"$s"/tmp/total_SNVs.list.all.final)
            Number_of_Polymorphic_sites_percentage=$(bc <<< "scale = 2; 100*$Number_of_Polymorphic_sites/$length") 
            total_number_of_mutation=$(grep "$gene_locus" results/13_snippy"$s"/tmp/core.vcf.2 | awk '{sum += $3} END {print sum}' ) ; if [ -z "$total_number_of_mutation" ] ; then total_number_of_mutation=$(echo "0") ; fi
            total_number_of_mutation_percentage=$(bc <<< "scale = 2; 100*$total_number_of_mutation/$length") ; 
            echo $gene_locus $gene_product $length $Number_of_Polymorphic_sites $Number_of_Polymorphic_sites_percentage $total_number_of_mutation $total_number_of_mutation_percentage >> results/13_snippy"$s"/tmp/gene_locus.Poly_sites.tab
        done
        ##---------------------------------------------------
        echo "running dN/dS"
        for F1 in $(cat $list); do
            (rm results/13_snippy"$s"/tmp/"$F1".dNdS.tmp)> /dev/null 2>&1 
            echo $F1 > results/13_snippy"$s"/"$F1"/tmp/"$F1".dNdS.tmp
            for gene in $(cat results/13_snippy"$s"/tmp/ref.gene.list); do
                echo $F1 $gene
                V1=$(grep "$gene" results/13_snippy"$s"/"$F1"/"$F1".13_snippy.SNPs-statistics.txt | awk '{print $6}') ; if [ -z "$V1" ] ; then V1=$(echo "0") ; fi #dN/dS
                #V2=$(grep "$gene" results/13_snippy"$s"/"$F1"/"$F1".13_snippy.SNPs-statistics.txt | awk '{print $3}') ; if [ -z "$V1" ] ; then V1=$(echo "0") ; fi #SNVs
                echo "$V1" >> results/13_snippy"$s"/tmp/"$F1".dNdS.tmp
                #echo "$V2" >> results/13_snippy"$s"/tmp/"$F1".SNV.tmp
            done
        done

        paste results/13_snippy"$s"/tmp/*.dNdS.tmp > results/13_snippy"$s"/tmp/all.dNdS.tmp
        #paste results/13_snippy"$s"/tmp/*.SNV.tmp > results/13_snippy"$s"/tmp/all.SNV.tmp

        awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' results/13_snippy"$s"/tmp/all.dNdS.tmp > results/13_snippy"$s"/tmp/all.dNdS.tmp.2
        awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' results/13_snippy"$s"/tmp/all.SNV.tmp > results/13_snippy"$s"/tmp/all.SNV.tmp.2

        paste results/13_snippy"$s"/tmp/gene_locus.Poly_sites.tab results/13_snippy"$s"/tmp/all.dNdS.tmp.2 results/13_snippy"$s"/tmp/all.SNV.tmp.2 > results/13_snippy"$s"/results.csv

        sed -i '1 i\gene_locus gene_product length Number_of_Polymorphic_sites $Number_of_Polymorphic_sites_percentage total_number_of_mutation total_number_of_mutation_percentage dNdS SNVs' results/13_snippy"$s"/results.csv
        sed -i "s/\t/ /g" results/13_snippy"$s"/results.csv

        ssconvert results/13_snippy"$s"/results.csv results/13_snippy"$s"/results.csv.xlsx
        ###############################################################################
        # calculate transition by tranversion ratio

        (java -jar /home/swapnil/pipeline/tools/SnpSift.jar tstv results/13_snippy"$s"/core.vcf > results/13_snippy"$s"/tmp/tstv.csv) > /dev/null 2>&1
        awk 'NR==3 || NR==6' results/13_snippy"$s"/tmp/tstv.csv > results/13_snippy"$s"/tmp/tstv.csv.tmp
        sed -i 's/,/\t/g' results/13_snippy"$s"/tmp/tstv.csv.tmp

        ##-------------------
        awk '
        { 
            for (i=1; i<=NF; i++)  {
                a[NR,i] = $i
            }
        }
        NF>p { p = NF }
        END {    
            for(j=1; j<=p; j++) {
                str=a[1,j]
                for(i=2; i<=NR; i++){
                    str=str" "a[i,j];
                }
                print str
            }
        }' results/13_snippy"$s"/tmp/tstv.csv.tmp > results/13_snippy"$s"/tmp/tstv.csv.transposed.tmp
        ##-------------------

        for F1 in $(cat $list); do
            A1=$(grep $F1 results/13_snippy"$s"/tmp/tstv.csv.transposed.tmp)
            echo $A1 > results/13_snippy"$s"/tmp/$F1.tstv.txt
        done

        rm results/13_snippy"$s"/tmp/tstv.csv.tmp
        rm results/13_snippy"$s"/tmp/tstv.csv.transposed.tmp
    fi
###############################################################################
echo "completed.... step-13 snippy -------------------------------------------"
###############################################################################
