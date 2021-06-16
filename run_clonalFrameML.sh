#!/bin/bash
###############################################################################
## ClonalFrameML doesnt work for roary output
## ClonalFrameML work for mugsy or mauve (or parsnp alignment)

if [ ! -d results/17_roary/ClonalFrameML ] ; then
mkdir results/17_roary/ClonalFrameML
fi

ClonalFrameML results/17_roary/extracted-snps.snp_sites.aln.fasttree.tree results/17_roary/roary_results/core_gene_alignment.aln results/17_roary/ClonalFrameML/ClonalFrameML_output
Rscript /home/swapnil/tools/ClonalFrameML-master/src/cfml_results.R results/17_roary/ClonalFrameML/ClonalFrameML_output

/home/swapnil/tools/maskrc-svg-master/maskrc-svg.py --aln results/17_roary/roary_results/core_gene_alignment.aln --out results/17_roary/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta results/17_roary/ClonalFrameML/ClonalFrameML_output

snp-sites -m -o results/17_roary/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta results/17_roary/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta

fasttree -nt -gtr < results/17_roary/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta > results/17_roary/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta.fasttree.tree

Core_genome_Rec_unfiltered_length=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_roary/roary_results/core_gene_alignment.aln | awk 'NR==1 {print $2}' )
Core_genome_Rec_unfiltered_SNVs=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_roary/extracted-snps.snp_sites.fasta | awk 'NR==1 {print $2}' )
Core_genome_Rec_filtered_SNVs=$(bioawk -c fastx '{ print $name, length($seq) }' <results/17_roary/ClonalFrameML/core_gene_alignment.aln.maskrc.fasta.snp_sites.fasta | awk 'NR==1 {print $2}' )

echo "$Core_genome_length $Core_genome_Rec_unfiltered_SNVs $Core_genome_Rec_filtered_SNVs" > results/17_roary/ClonalFrameML/stat.tab
