###############################################################################
    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate myenv
###############################################################################
mkdir results/17_roary_all/coinfinder_results

coinfinder -i results/17_roary_all/roary_results/gene_presence_absence.csv -I -p results/17_roary_all/extracted-snps.snp_sites.aln.fasttree.tree -o results/17_roary_all/coinfinder_results --associate -x 8
coinfinder -i results/17_roary_all/roary_results/gene_presence_absence.csv -I -p results/17_roary_all/extracted-snps.snp_sites.aln.fasttree.tree -o results/17_roary_all/coinfinder_results --dissociate -x 8 

###############################################################################