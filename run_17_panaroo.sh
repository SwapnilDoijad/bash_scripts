#!/bin/bash
###############################################################################
#9 WGS-alignment-mugsy
# add -z option to keep the intermediate files
###############################################################################
echo "started.... step-17 panaroo ----------------------------------------------"
###############################################################################
echo "need to modify the script"
## initial inout and file preparation
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | sed 's/list\.//g' | sed 's/\.txt//g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "want to suffix files? type "y" else press enter to continue"
    read answer
    if [ "$answer" = "y" ]; then
    s="_$l"
    fi

    echo "Wish to run ClonalFrameML?"
    read CFML_ans

    (mkdir results/17_panaroo"$s") > /dev/null 2>&1
    (mkdir results/17_panaroo"$s"/gff) > /dev/null 2>&1

    for F1 in $(cat $list);do
    cp results/08_annotation/raw_files/$F1/$F1.gff results/17_panaroo"$s"/gff/
    done
###############################################################################
## run panaroo

    #source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    #conda activate

    ##if conda is not working, then use panaroo installed by pip3 (pip3 install git+https://github.com/gtonkinhill/panaroo)
    (mkdir results/17_panaroo"$s"/panaroo_results) > /dev/null 2>&1

    panaroo \
    -i results/17_panaroo"$s"/gff/*.gff \
    -o results/17_panaroo"$s"/panaroo_results \
    -t 8 \
    -c 0.70 \
    --len_dif_percent 0.7 \
    --quiet \
    --clean-mode strict \
    --merge_paralogs \
    -a core

    #conda deactivate

exit
###############################################################################

    # snp-sites -m -p -o results/17_panaroo"$s"/extracted-snps  results/17_panaroo"$s"/panaroo_results/core_gene_alignment.aln

    # (fasttree -nt -gtr < results/17_panaroo"$s"/extracted-snps.snp_sites.aln > results/17_panaroo"$s"/extracted-snps.snp_sites.aln.fasttree.tree)> /dev/null 2>&1

    # (fastme -i results/17_panaroo"$s"/extracted-snps.phylip -o results/17_panaroo"$s"/tree.fastme.nwk -b 100 -n -d)> /dev/null 2>&1

    # (raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s results/17_panaroo"$s"/extracted-snps.phylip -n T20 )> /dev/null 2>&1

    mv RAxML_bipartitions.T20 results/17_panaroo"$s"/

    rm *.T20
    rm *.T19

    rm *.phy_fastme_stat.txt
    rm *.phy_fastme_boot.txt

    Core_genome_Rec_unfiltered_length=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_panaroo"$s"/panaroo_results/core_gene_alignment.aln | awk 'NR==1 {print $2}' )
    Core_genome_Rec_unfiltered_SNVs=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_panaroo"$s"/extracted-snps.snp_sites.fasta | awk 'NR==1 {print $2}' )

    echo "$Core_genome_length $Core_genome_Rec_unfiltered_SNVs" > results/17_panaroo"$s"/ClonalFrameML/stat.tab
    rm -rf results/17_panaroo"$s"/gff
    echo "completed.... step-16 panaroo --------------------------------------------"

###############################################################################
##Run ClonanFrameML
    if [ "$CFML_ans" = "y" ] ; then
        echo "running clonalFrameML"
        (mkdir results/17_panaroo"$s"/ClonalFrameML) > /dev/null 2>&1
        ClonalFrameML results/17_panaroo"$s"/extracted-snps.snp_sites.aln.fasttree.tree results/17_panaroo"$s"/panaroo_results/core_gene_alignment.aln results/17_panaroo"$s"/ClonalFrameML/ClonalFrameML_output
        Rscript /home/swapnil/tools/ClonalFrameML-master/src/cfml_results.R results/17_panaroo"$s"/ClonalFrameML/ClonalFrameML_output
        /home/swapnil/tools/maskrc-svg-master/maskrc-svg.py --aln results/17_panaroo"$s"/panaroo_results/core_gene_alignment.aln --out results/17_panaroo"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta results/17_panaroo"$s"/ClonalFrameML/ClonalFrameML_output
        snp-sites -m -o results/17_panaroo"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta results/17_panaroo"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta
        fasttree -nt -gtr < results/17_panaroo"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta > results/17_panaroo"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta.fasttree.tree
        Core_genome_Rec_filtered_SNVs=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_panaroo"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta | awk 'NR==1 {print $2}' )
        echo "$Core_genome_length $Core_genome_Rec_unfiltered_SNVs $Core_genome_Rec_filtered_SNVs" > results/17_panaroo"$s"/ClonalFrameML/stat.tab
    fi
###############################################################################

