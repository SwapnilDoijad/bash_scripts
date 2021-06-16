###############################################################################
echo "started.... funannotate ------------------------------------------------"
###############################################################################
#------------------------------------------------------------------------------
## installation step

    ## check that all modules are installed
    #funannotate check --show-versions

    ##download/setup databases to a writable/readable location
    # funannotate setup -d $HOME/funannotate_db

    ##set ENV variable for $FUNANNOTATE_DB
    #echo "export FUNANNOTATE_DB=$HOME/funannotate_db" > /conda/installation/path/envs/funannotate/etc/conda/activate.d/funannotate.sh
    #echo "unset FUNANNOTATE_DB" > /conda/installation/path/envs/funannotate/etc/conda/deactivate.d/funannotate.sh

    ##run tests -- requires internet connection to download data
    #funannotate test -t all --cpus 16

#------------------------------------------------------------------------------
###############################################################################
## funannotate
source /home/swapnil/miniconda3/etc/profile.d/conda.sh
conda activate funannotate

echo "provide list file (for e.g. all)"
echo "-------------------------------------------------------------------------"
ls list.*.txt | sed 's/ /\n/g'
echo "-------------------------------------------------------------------------"
read l
list=$(echo "list.$l.txt")

echo "provide species name"
read spp

echo "do you wish funannotate_masking? type y for yes"
read funannotate_masking

echo "do you wish funannotate_masked contigs to join? type y for yes"
read contigs_join

annot_dir=results/08_annotation
assembly_dir=results/04_assembly

(mkdir results)> /dev/null 2>&1
(mkdir $annot_dir)> /dev/null 2>&1

###############################################################################
for F1 in $(cat $list); do
    (mkdir $assembly_dir/$F1)> /dev/null 2>&1
    (mkdir $assembly_dir/$F1/funannotate)> /dev/null 2>&1
    work_dir=$assembly_dir/$F1/funannotate

    ## funannotate masking ----------------------------------------------------
    ## funannotate cleaning, sorting and masking 
    if [ "$funannotate_masking" == "y" ]; then
        if [ -f $assembly_dir/all_fasta/$F1.tantan_masked.fasta ]; then
            echo "funannotate cleaning, sorting and masking  for $F1 is already finished"
            else
            echo "running funannotate cleaning, sorting and masking for $F1"
            funannotate clean -i $assembly_dir/$F1/$F1.fasta -o $work_dir/$F1.cleaned.fasta
            funannotate sort -i $work_dir/$F1.cleaned.fasta -o $work_dir/$F1.cleaned.sorted.fasta
            funannotate mask -i $work_dir/$F1.cleaned.sorted.fasta -o $work_dir/$F1.cleaned.sorted.masked.fasta
            cp $work_dir/$F1.cleaned.sorted.masked.fasta $assembly_dir/all_fasta/$F1.tantan_masked.fasta
        fi
    fi

    ## funannotate contig joining step
    if [ "$contigs_join" == "y" ]; then
        if [ -f $assembly_dir/all_fasta/$F1.tantan_masked.joined.fasta ]; then 
            echo "funannotate contig joining step for $F1 already finished"
            else
            echo "running funannotate contig joining step for $F1"
            cp $work_dir/$F1.cleaned.sorted.masked.fasta $assembly_dir/$F1/$F1.cleaned.sorted.masked.joined.fasta.tmp
            sed -i 's/>.*/NNNNNNNNNN/g' $assembly_dir/$F1/$F1.cleaned.sorted.masked.joined.fasta.tmp
            sed -i "1i "'>'$F1"" $assembly_dir/$F1/$F1.cleaned.sorted.masked.joined.fasta.tmp
            fasta_formatter -i $assembly_dir/$F1/$F1.cleaned.sorted.masked.joined.fasta.tmp -w 60 -o $assembly_dir/$F1/$F1.cleaned.sorted.masked.joined.fasta
            rm $assembly_dir/$F1/$F1.cleaned.sorted.masked.joined.fasta.tmp
            cp $assembly_dir/$F1/$F1.cleaned.sorted.masked.joined.fasta $assembly_dir/all_fasta/$F1.tantan_masked.joined.fasta
        fi
    fi

    ## funannotate prediction step for non-joined fasta
    if [ -f results/04_assembly/all_fasta/$F1.tantan_masked.fasta ] ; then
        if [ -d results/08_annotation/$F1/funannotate_tantan-masked_fasta ]; then
            echo "Funnannotate CDS prediction for $F1 tantan_masked.fasta (non-joined fasta) already finished"
            else
            echo "running CDS prediction Funnannotate for $F1 tantan_masked.fasta (non-joined fasta)"
            (mkdir results/08_annotation/$F1/funannotate_tantan-masked_fasta )> /dev/null 2>&1
            funannotate predict -i results/04_assembly/all_fasta/$F1.tantan_masked.fasta -o results/08_annotation/$F1/funannotate_tantan-masked_fasta -s "$spp"
        fi
    fi

    ## funannotate prediction step for joined fasta
    if [ "$contigs_join" == "y" ]; then
        if [ -d results/08_annotation/$F1/funannotate_tantan-masked_joined_fasta ]; then
            echo "funannotate prediction step for $F1 joined fasta already finished"
            else
            echo "running funannotate prediction step for $F1 joined fasta"
            if [ -f results/04_assembly/all_fasta/$F1.tantan_masked.joined.fasta ] ; then
                echo "running Funnannotate for $F1 tantan_masked.joined.fasta"
                (mkdir results/08_annotation/$F1/funannotate_tantan-masked_joined_fasta )> /dev/null 2>&1
                funannotate predict -i results/04_assembly/all_fasta/$F1.tantan_masked.joined.fasta -o results/08_annotation/$F1/funannotate_tantan-masked_joined_fasta -s "$spp"
            fi
        fi
    fi

    ## RepeatScout masked fasta -----------------------------------------------    
    ## CDS prediction step for non-joined fasta
    if [ -f results/04_assembly/all_fasta/$F1.repeatscout_masked.fasta ] ; then
        if [ -d results/08_annotation/$F1/funannotate_repeatscout-masked_fasta ]; then
            echo "RepeatScout masked fasta: CDS predictoin step for $F1 already finished (non-joined)"
            else
            echo "RepeatScout masked fasta: running CDS predictoin step for $F1 (non-joined)"
            (mkdir results/08_annotation/$F1/funannotate_repeatscout-masked_fasta )> /dev/null 2>&1
            funannotate predict -i results/04_assembly/all_fasta/$F1.repeatscout_masked.fasta -o results/08_annotation/$F1/funannotate_repeatscout-masked_fasta -s "$spp"
        fi
    fi
    ## CDS prediction step for joined fasta
    if [ "$contigs_join" == "y" ]; then
        if [ -f results/04_assembly/all_fasta/$F1.repeatscout_masked.joined.fasta ] ; then
            if [ -d results/08_annotation/$F1/funannotate_repeatscout-masked_joined_fasta ]; then
                echo "RepeatScout masked fasta: running Funnannotate for $F1 repeatscout_masked.joined.fasta"
                else
                echo "RepeatScout masked fasta: running Funnannotate for $F1 repeatscout_masked.joined.fasta"
                (mkdir results/08_annotation/$F1/funannotate_repeatscout-masked_joined_fasta )> /dev/null 2>&1
                funannotate predict -i results/04_assembly/all_fasta/$F1.repeatscout_masked.joined.fasta -o results/08_annotation/$F1/funannotate_repeatscout-masked_joined_fasta -s "$spp"
            fi
        fi
    fi

done
conda deactivate
###############################################################################
echo "finished... funannotate ------------------------------------------------"
###############################################################################
## tantan-masked fasta

## 
    ##eggnog mapping
    conda activate funannotate ##required for diamond version conflict
    /home/swapnil/tools/eggnog-mapper-master/emapper.py \
    -i results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results/Alternaria_alternata.proteins.fa \
    --output Alternaria_alternata.proteins.fa.diamond.eggnog \
    --output_dir results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results/eggnog/ \
    -m diamond --dmnd_db /media/swapnil/databases/bioinfoDBs/eggnog_database/eggnog_proteins.dmnd \
    --data_dir /media/swapnil/databases/bioinfoDBs/eggnog_database/ \
    --cpu 8 
    conda activate funannotate

    ## interproscan
    /home/swapnil/tools/interproscan-5.19-58.0/interproscan.sh \
    -i results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results_eggnog/Alternaria_alternata.proteins.fa \
    -b Alternaria_alternata.proteins.fa.interproscan \
    -f xml \
    -goterms \
    -iprlookup \
    -pa 

    mv Alternaria_alternata.proteins.fa.interproscan.xml results/08_annotation/alternaria/funannotate_tantan-masked_fasta/interproscan
    ###############################################################################
    ## antismash
    conda activate antismash
        antismash results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results_eggnog/Alternaria_alternata.gbk --taxon fungi 
        #mv results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results_eggnog/$F1 results/35_antiSMASH/results/
    conda deactivate antismash
    ###############################################################################
    ## phobius
    ## not imp phobius results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results_eggnog/Alternaria_alternata.proteins.fa

    ###############################################################################
    ## final annotation
    funannotate annotate \
    -i results/08_annotation/alternaria/funannotate_tantan-masked_fasta/predict_results \
    -o funannotate_tantan-masked_fasta_functional_annotation \
    --fasta results/04_assembly/alternaria/funannotate/alternaria.cleaned.sorted.tantan-masked.fasta \
    -s "Alternaria alternata" \
    --eggnog results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results/eggnog/Alternaria_alternata.proteins.fa.diamond.eggnog.emapper.annotations \
    --antismash results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results/antismash/scaffold_1.final.gbk \
    --iprscan results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results/interproscan/Alternaria_alternata.proteins.fa.interproscan.xml









 ###############################################################################
    ## RepeatScout-masked fasta

## 
    ##eggnog mapping

    conda activate funannotate ##required for diamond version conflict
    /home/swapnil/tools/eggnog-mapper-master/emapper.py \
    -i results/08_annotation/alternaria/funannotate_repeatscout-masked_fasta/predict_results/Alternaria_alternata.proteins.fa \
    --output Alternaria_alternata.proteins.fa.diamond.eggnog \
    --output_dir results/08_annotation/alternaria/funannotate_repeatscout-masked_fasta/predict_results/eggnog/ \
    -m diamond --dmnd_db /media/swapnil/databases/bioinfoDBs/eggnog_database/eggnog_proteins.dmnd \
    --data_dir /media/swapnil/databases/bioinfoDBs/eggnog_database/ \
    --cpu 8
    conda deactivate


    ## interproscan
    /home/swapnil/tools/interproscan-5.19-58.0/interproscan.sh \
    -i results/08_annotation/alternaria/funannotate_repeatscout-masked_fasta/predict_results/Alternaria_alternata.proteins.fa \
    -b Alternaria_alternata.proteins.fa.interproscan \
    -f xml \
    -goterms \
    -iprlookup \
    -pa 

    mv Alternaria_alternata.proteins.fa.interproscan.xml results/08_annotation/alternaria/funannotate_tantan-masked_fasta/interproscan
    ###############################################################################
    ## antismash
    conda activate antismash
        antismash results/08_annotation/alternaria/funannotate_repeatscout-masked_fasta/predict_results/Alternaria_alternata.gbk --taxon fungi
        #mv results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results_eggnog/$F1 results/35_antiSMASH/results/
    conda deactivate antismash
    ###############################################################################
    ## phobius
    ## not imp phobius results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results_eggnog/Alternaria_alternata.proteins.fa

    ###############################################################################
    ## final annotation
    funannotate annotate \
    -i results/08_annotation/alternaria/funannotate_tantan-masked_fasta/predict_results \
    -o funannotate_tantan-masked_fasta_functional_annotation \
    --fasta results/04_assembly/alternaria/funannotate/alternaria.cleaned.sorted.tantan-masked.fasta \
    -s "Alternaria alternata" \
    --eggnog results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results/eggnog/Alternaria_alternata.proteins.fa.diamond.eggnog.emapper.annotations \
    --antismash results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results/antismash/scaffold_1.final.gbk \
    --iprscan results/08_annotation/alternaria/funannotate_tantan-masked_fasta/annotate_results/interproscan/Alternaria_alternata.proteins.fa.interproscan.xml