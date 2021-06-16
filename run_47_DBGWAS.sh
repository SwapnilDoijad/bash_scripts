###############################################################################
# 47 DBGWAS
## If strain numbers are too less, program fails !!!!
###############################################################################
## Preliminary input and directories creation
    ncdb=/home/swapnil/tools/DBGWAS-0.5.4-Linux-precompiled/Resistance_DB_for_DBGWAS.fasta
    pcdb=/home/swapnil/tools/DBGWAS-0.5.4-Linux-precompiled/uniprot_sprot_bacteria_for_DBGWAS.fasta

#    echo -e "place data/traits.dbgwas.TRAIT.txt(unrelated strains should be labeled as NA) file "

    echo "provide trait name (for e.g. gentamicin)"
    echo "-------------------------------------------------------------------------"
    ls data/traits.*.txt | sed 's/data\/traits\.dbgwas\.//g' | sed 's/\.txt//g' | sed 's/ /\n/g' 
    echo "-------------------------------------------------------------------------"
    read l
    list=$(awk -F',' '{print $1}' data/traits.dbgwas.$l.txt)
    s=$(echo "_$l")

    (mkdir results) > /dev/null 2>&1
    (mkdir results/47_DBGWAS"$s") > /dev/null 2>&1
    (mkdir results/47_DBGWAS"$s"/fasta) > /dev/null 2>&1
    (mkdir results/47_DBGWAS"$s"/tmp) > /dev/null 2>&1

#    echo "MIC? (for e.g 4) applicable for MICs, else just press enter"
#    read MIC
###############################################################################
    echo "ID Phenotype Path" > results/47_DBGWAS"$s"/pheno.2.txt
    for F1 in $list ; do
        V1=$(grep "$F1" data/traits.dbgwas.$l.txt | sed 's/,/\t/g')
        V2=$(echo "results/47_DBGWAS"$s"/fasta/$F1.fasta")
        echo $V1 $V2 >> results/47_DBGWAS"$s"/pheno.2.txt
        cp results/04_assembly/all_fasta/$F1.fasta results/47_DBGWAS"$s"/fasta/
    done
    sed -i 's/ /\t/g' results/47_DBGWAS"$s"/pheno.2.txt
    mash sketch -s 10000 -o results/47_DBGWAS"$s"/tmp/mash_sketch results/47_DBGWAS"$s"/fasta/*.fasta
    mash dist results/47_DBGWAS"$s"/tmp/mash_sketch.msh results/47_DBGWAS"$s"/tmp/mash_sketch.msh | square_mash > results/47_DBGWAS"$s"/tmp/mash.tsv
    mv  results/47_DBGWAS"$s"/tmp/mash.tsv results/mash.tsv
    Rscript /home/swapnil/pipeline/tools/plot_and_tree.47_DBGWAS.r 
    mv results/tree.nwk results/47_DBGWAS"$s"/tmp/tree.nwk
    mv results/mash.tsv  results/47_DBGWAS"$s"/tmp/mash.tsv
    
    rm -rf results/47_DBGWAS"$s"/output
    
    cp data/traits.dbgwas.$l.txt results/47_DBGWAS"$s"/traits.dbgwas.$l.txt
    echo "running DBGWAS for $l"
    DBGWAS -strains results/47_DBGWAS"$s"/pheno.2.txt -newick results/47_DBGWAS"$s"/tmp/tree.nwk -nc-db $ncdb -pt-db $pcdb -output results/47_DBGWAS"$s"/output -nb-cores 8 
    
    rm -rf results/47_DBGWAS"$s"/fasta/ ## delete to save space and avoid redundancy
    rm -rf results/47_DBGWAS"$s"/output/step1 ## step1 folder very big in size, so delete
    ########################################### other commands
    ## if phenotype threshold is mentioned
    ##DBGWAS -strains results/47_DBGWAS"$s"/pheno.2.txt -newick results/47_DBGWAS"$s"/tmp/tree.nwk -nc-db $ncdb -pt-db $pcdb -output results/47_DBGWAS"$s"/output -phenoThreshold $MIC -nb-cores 8
    ## with -Rscript
    #DBGWAS -strains results/47_DBGWAS"$s"/pheno.2.txt -newick results/47_DBGWAS"$s"/tmp/tree.nwk -nc-db $ncdb -pt-db $pcdb -Rscript-path /home/swapnil/Downloads/R-3.4.0/bin/Rscript -output results/47_DBGWAS"$s"/output -nb-cores 8 
    ## normal command-line
    ## docker version
    # sudo docker run --rm -v /home/swapnil/p_Enterobacter_NCBI/results/47_DBGWAS_test:/47_DBGWAS_test leandroishilima/dbgwas:0.5.4 -strains 47_DBGWAS"$s"/pheno.2.txt -newick 47_DBGWAS"$s"/tmp/tree.nwk -nc-db $ncdb -pt-db $pcdb -output 47_DBGWAS"$s"/output -nb-cores 8
    # this one could run succesfully 
    # sudo docker run --rm -v /home/swapnil/p_Enterobacter_NCBI/results/47_DBGWAS_test:/47_DBGWAS_test leandroishilima/dbgwas:0.5.4 -strains 47_DBGWAS_test/pheno.2.txt -newick 47_DBGWAS_test/tmp/tree.nwk -output 47_DBGWAS_test/output -nb-cores 8
    
###############################################################################