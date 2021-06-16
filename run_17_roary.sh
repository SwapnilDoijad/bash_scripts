#!/bin/bash
###############################################################################
#9 WGS-alignment-mugsy
# add -z option to keep the intermediate files
###############################################################################
echo "started.... step-17 roary ----------------------------------------------"
###############################################################################
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

    (mkdir results/17_roary"$s") > /dev/null 2>&1
    (mkdir results/17_roary"$s"/gff) > /dev/null 2>&1

    for F1 in $(cat $list);do
    cp results/08_annotation/raw_files/$F1/$F1.gff results/17_roary"$s"/gff/
    done
###############################################################################
## run roary

    #(roary -p 8 -f results/17_roary"$s"/roary_results -e -n -z -i 70 -r results/17_roary"$s"/gff/*.gff ) > /dev/null 2>&1
    ( roary -p 8 -f results/17_roary"$s"/roary_results -e -n results/17_roary"$s"/gff/*.gff ) > /dev/null 2>&1
    #roary -p 8 -f results/17_roary"$s"/roary_results -g 75000 -e -n -i 70 -r results/17_roary"$s"/gff/*.gff 
    #query_pan_genome -g results/17_roary"$s"/roary_results/clustered_proteins -a difference --input_set_one results/17_roary"$s"/gff/EB-247_WGS.gff --input_set_two results/17_roary"$s"/gff/Survcare220.gff

    snp-sites -m -p -o results/17_roary"$s"/extracted-snps  results/17_roary"$s"/roary_results/core_gene_alignment.aln

    ## run FastTree (Note: sometime fasttree out put is not taken by clonalframeML)
    (fasttree -nt -gtr < results/17_roary"$s"/extracted-snps.snp_sites.aln > results/17_roary"$s"/extracted-snps.snp_sites.aln.fasttree.tree)> /dev/null 2>&1

    #(fastme -i results/17_roary"$s"/extracted-snps.phylip -o results/17_roary"$s"/tree.fastme.nwk -b 100 -n -d)> /dev/null 2>&1

    #(raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s results/17_roary"$s"/extracted-snps.phylip -n T20)> /dev/null 2>&1

    mv RAxML_bipartitions.T20 results/17_roary"$s"/

    rm *.T20
    rm *.T19

    rm *.phy_fastme_stat.txt
    rm *.phy_fastme_boot.txt

    Core_genome_Rec_unfiltered_length=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_roary"$s"/roary_results/core_gene_alignment.aln | awk 'NR==1 {print $2}' )
    Core_genome_Rec_unfiltered_SNVs=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_roary"$s"/extracted-snps.snp_sites.fasta | awk 'NR==1 {print $2}' )

    echo "$Core_genome_length $Core_genome_Rec_unfiltered_SNVs" > results/17_roary"$s"/ClonalFrameML/stat.tab
    rm -rf results/17_roary"$s"/gff
    echo "completed.... step-16 roary --------------------------------------------"
###############################################################################
## piggy
#    echo "running piggy"
#    /home/swapnil/tools/piggy/bin/piggy -i results/17_roary"$s"/gff/ -o results/17_roary"$s"/piggy -r results/17_roary"$s"/roary_results/ -t 8 -R
###############################################################################
## roary2fripan.py
#    /home/swapnil/tools/roary/roary2fripan-master/roary2fripan.py --input results/17_roary"$s"/roary_results/gene_presence_absence.csv test123
###############################################################################
##Run ClonanFrameML
    ## (Note: sometime fasttree out put is not taken by clonalframeML)
    if [ "$CFML_ans" = "y" ] ; then
        echo "running clonalFrameML"
        (mkdir results/17_roary"$s"/ClonalFrameML) > /dev/null 2>&1
        ClonalFrameML results/17_roary"$s"/extracted-snps.snp_sites.aln.fasttree.tree results/17_roary"$s"/roary_results/core_gene_alignment.aln results/17_roary"$s"/ClonalFrameML/ClonalFrameML_output
        Rscript /home/swapnil/tools/ClonalFrameML-master/src/cfml_results.R results/17_roary"$s"/ClonalFrameML/ClonalFrameML_output
        /home/swapnil/tools/maskrc-svg-master/maskrc-svg.py --aln results/17_roary"$s"/roary_results/core_gene_alignment.aln --out results/17_roary"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta results/17_roary"$s"/ClonalFrameML/ClonalFrameML_output
        snp-sites -m -o results/17_roary"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta results/17_roary"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta
        fasttree -nt -gtr < results/17_roary"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta > results/17_roary"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta.fasttree.tree
        Core_genome_Rec_filtered_SNVs=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_roary"$s"/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta | awk 'NR==1 {print $2}' )
        echo "$Core_genome_length $Core_genome_Rec_unfiltered_SNVs $Core_genome_Rec_filtered_SNVs" > results/17_roary"$s"/ClonalFrameML/stat.tab
    fi
###############################################################################
## run fastGEAR for recombination
    ## Note: provide full path, otherwise sometime doesenot work
    echo "running fastGEAR"
    mkdir results/52_fastGEAR_ancestral"$s"
    bash /home/swapnil/tools/fastGEARpackageLinux64bit/run_fastGEAR.sh /usr/local/MATLAB/MATLAB_Runtime/v901 results/17_roary"$s"/roary_results/core_gene_alignment.aln results/52_fastGEAR_ancestral"$s"/out.mat /home/swapnil/tools/fastGEARpackageLinux64bit/fG_input_specs.txt
###############################################################################