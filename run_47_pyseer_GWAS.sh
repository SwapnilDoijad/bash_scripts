#!/bin/bash

## use bcftools to extact information https://samtools.github.io/bcftools/howtos/query.html
## for e.g 
## bcftools query -f '%POS,%INFO/TYPE\n' merged.vcf
## bcftools query -f '%POS,%REF,%ALT,%INFO/TYPE,%INFO/ANN\n' merged.vcf | sed 's/|/\t/g' > test.csv
## bcftools query -f '%POS,%REF,%ALT,%INFO/TYPE,%INFO/ANN\n' merged.vcf | sed 's/|/\t/g' > merged.vcf.info.csv

## conda pyseer not working, use system default installed pyseer
###############################################################################
# 47 GWAS
###############################################################################
## Preliminary input and directories creation
    (mkdir results) > /dev/null 2>&1
    (mkdir results/47_GWAS_pyseer) > /dev/null 2>&1
    (mkdir results/47_GWAS_pyseer/fasta) > /dev/null 2>&1
    (mkdir results/47_GWAS_pyseer/tmp) > /dev/null 2>&1
    (mkdir results/47_GWAS_pyseer/snippy) > /dev/null 2>&1
    echo -e "place GWAS_pheno.txt(first line is heading, unrelated strains should be labeled as NA) file in results/GWAS folder"
    echo "provide list file (for e.g. all)"
    ls list.*.txt | sed 's/ /\n/g'
    read l
    list=$(echo "list.$l.txt")
    s=$(echo "_$l")
    strain=$(cat $list | wc -l)
    strains_10percentage=$(( $strain * 10/100 ))
    strains=$(( $strain - $strains_10percentage ))
    echo "provide reference with path (for e.g. results/00_ref/ref.ESBL3012.gbk)"
    ls results/00_ref/*.gbk | sed 's/ /\n/g'
    read ref_path
    ref=$(echo $ref_path | awk -F'/' '{print $NF}' | sed 's/.gbk//g' | sed 's/ref.//g')
    echo "provide path for gene_presence_absence.Rtab (for e.g. results/17_roary*/roary_results/gene_presence_absence.Rtab)"
    find results/17_roary* -name "gene_presence_absence.Rtab"
    read gpa
    echo "provide snippy path (for e.g. results/13_snippy_Ex_84)"
    ls results/ | grep 13_snippy
    read snippy_path_tmp
    snippy_path=$(echo "results/$snippy_path_tmp")
###############################################################################
## copy fasta, calculate mash and run scree_plot
    echo "running mash and scree_plot"
    for F1 in $(cat $list); do
        cp results/04_assembly/all_fasta/$F1.fasta results/47_GWAS_pyseer/fasta/
    done
    mash sketch -s 10000 -o results/47_GWAS_pyseer/tmp/mash_sketch results/47_GWAS_pyseer/fasta/*.fasta
    mash dist results/47_GWAS_pyseer/tmp/mash_sketch.msh results/47_GWAS_pyseer/tmp/mash_sketch.msh | square_mash > results/47_GWAS_pyseer/tmp/mash.tsv
    Rscript /home/swapnil/pipeline/tools/plot_and_tree.47_DBGWAS.r
    ## python scripts/phylogeny_distance.py core_genome_aln.tree > phylogeny_dists
    (scree_plot_pyseer results/47_GWAS_pyseer/tmp/mash.tsv) > /dev/null 2>&1
    mv scree_plot.png results/47_GWAS_pyseer/tmp/
    echo "check the scree_plot and provide the max-dimensions values"
    read md
###############################################################################
## run the COGs analysis:
    echo "running GWAS for COGs"    
    pyseer --phenotypes data/GWAS_pheno.txt --pres $gpa --distances results/47_GWAS_pyseer/tmp/mash.tsv --save-m results/47_GWAS_pyseer/tmp/mash_mds --max-dimensions $md > results/47_GWAS_pyseer/tmp/phenotype_COGs.txt
    sort -t$'\t' -k4,4 results/47_GWAS_pyseer/phenotype_COGs.sorted.csv
###############################################################################
## snp based analysis 
    echo "running GWAS for SNPs"  
    for F1 in $(cat $list); do
    cp $snippy_path/$F1/$F1.vcf.gz results/47_GWAS_pyseer/snippy/
    cp $snippy_path/$F1/$F1.vcf.gz.csi results/47_GWAS_pyseer/snippy/
    done
    bcftools merge -m none -0 -O z results/47_GWAS_pyseer/snippy/*.vcf.gz > results/47_GWAS_pyseer/tmp/merged.vcf.gz
    pyseer --phenotypes data/GWAS_pheno.txt --vcf results/47_GWAS_pyseer/tmp/merged.vcf.gz --load-m results/47_GWAS_pyseer/tmp/mash_mds.pkl --min-af 0.02 --max-af 0.98 > results/47_GWAS_pyseer/phenotype_SNPs.txt
    #pyseer --phenotypes data/GWAS_pheno.txt --vcf results/47_GWAS_pyseer/tmp/merged.vcf.gz --load-m results/47_GWAS_pyseer/tmp/mash_mds.pkl --lineage --print-samples > results/47_GWAS_pyseer/tmp/phenotype_SNPs.txt
    ## chek below if the reference fasta has a specical chracter
    cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' results/47_GWAS_pyseer/tmp/phenotype_SNPs.txt | cut -d "_" -f 2) <(sed '1d' results/47_GWAS_pyseer/tmp/phenotype_SNPs.txt | cut -f 4) | awk '{p = -log($2)/log(10); print "chr",".",$1,p,p,"0"}' ) | tr ' ' '\t' > results/47_GWAS_pyseer/tmp/phenotype_snps.plot
###############################################################################
## kmer based analysis
    echo "runnign through kmer analysis"
    (rm results/47_GWAS_pyseer/fsm_file_list.txt) > /dev/null 2>&1
    for F1 in $(cat $list); do
    echo "$F1 results/47_GWAS_pyseer/fasta/$F1.fasta" | sed 's/ /\t/g' >> results/47_GWAS_pyseer/fsm_file_list.txt
    done
    /home/swapnil/tools/fsm-lite-master/fsm-lite -l results/47_GWAS_pyseer/fsm_file_list.txt -s $strains_10percentage -S $strains -v -t fsm_kmers | gzip -c - > results/47_GWAS_pyseer/fsm_kmers.txt.gz
    python3.6 /home/swapnil/pipeline/tools/phylogeny_distance.py --lmm results/47_GWAS_pyseer/tmp/tree.nwk > results/47_GWAS_pyseer/tmp/phylogeny_K.tsv
    pyseer --lmm --phenotypes data/GWAS_pheno.txt --kmers results/47_GWAS_pyseer/fsm_kmers.txt.gz --similarity results/47_GWAS_pyseer/tmp/phylogeny_K.tsv --output-patterns kmer_patterns.txt --cpu 16 > results/47_GWAS_pyseer/tmp/phenotype_kmers.txt
    python3.6 /home/swapnil/pipeline/tools/count_patterns.py results/47_GWAS_pyseer/tmp/phenotype_kmers.txt > results/47_GWAS_pyseer/tmp/phenotype_kmers.patterns.txt
    python3.6 /home/swapnil/pipeline/tools/qq_plot.py results/47_GWAS_pyseer/tmp/phenotype_kmers.txt
    mv qq_plot.png results/47_GWAS_pyseer/tmp/
    cat <(head -1 results/47_GWAS_pyseer/tmp/phenotype_kmers.txt) <(awk '$4<'1E-04' {print $0}' results/47_GWAS_pyseer/tmp/phenotype_kmers.txt) > results/47_GWAS_pyseer/tmp/significant_kmers.txt
    phandango_mapper results/47_GWAS_pyseer/tmp/significant_kmers.txt results/00_ref/$ref.fna results/47_GWAS_pyseer/tmp/kmers.plot

    ## annotating kmer
    annotate_hits_pyseer results/47_GWAS_pyseer/tmp/significant_kmers.txt results/47_GWAS_pyseer/tmp/references.txt results/47_GWAS_pyseer/tmp/annotated_kmers.txt
    python3.6 /home/swapnil/pipeline/tools/summarise_annotations.py results/47_GWAS_pyseer/tmp/annotated_kmers.txt > results/47_GWAS_pyseer/tmp/gene_hits.txt
    cp /home/swapnil/pyseer_Eb_2/results/GWAS/tmp/gene_hits.txt 
    Rscript  /home/swapnil/pyseer_Eb_2/results/GWAS/tmp/pyseer.ggplot.r
    mv Rplot.pdf results/47_GWAS_pyseer/
###############################################################################
# GWAS end
###############################################################################