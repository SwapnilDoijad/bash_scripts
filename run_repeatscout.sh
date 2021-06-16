###############################################################################
echo "started.... repeatscout and RepeatMasker --------------------------------"
###############################################################################
echo "provide list file (for e.g. all)"
echo "-------------------------------------------------------------------------"
ls list.*.txt | sed 's/ /\n/g'
echo "-------------------------------------------------------------------------"
read l
list=$(echo "list.$l.txt")

repeatscout_path=/home/swapnil/tools/RepeatScout-1.0.5/RepeatScout-1
repeatmasker_path=/home/swapnil/tools/RepeatMasker-4.1.0/RepeatMasker
assembly_dir=results/04_assembly

(mkdir results)> /dev/null 2>&1
(mkdir $assembly_dir)> /dev/null 2>&1

###############################################################################
for F1 in $(cat $list); do
    (mkdir $assembly_dir/$F1)> /dev/null 2>&1
    (mkdir $assembly_dir/$F1/repeatscout)> /dev/null 2>&1
    (mkdir $assembly_dir/$F1/repeatmasker)> /dev/null 2>&1
    cp $assembly_dir/$F1/$F1.fasta $assembly_dir/$F1/$F1.original.fasta #to save the original copy

    awk '/>/{sub(/>.*/,">contig_"++i)}1' $assembly_dir/$F1/$F1.fasta > $assembly_dir/$F1/$F1.fasta.tmp 
    rm $assembly_dir/$F1/$F1.fasta
    fasta_formatter -i $assembly_dir/$F1/$F1.fasta.tmp -w 60 -o $assembly_dir/$F1/$F1.fasta
    rm $assembly_dir/$F1/$F1.fasta.tmp

    echo "running step-1:build_lmer_table for $F1"
    $repeatscout_path/build_lmer_table -sequence $assembly_dir/$F1/$F1.fasta -freq $assembly_dir/$F1/repeatscout/$F1.output_lmer.frequency
    echo "running step-2:RepeatScout for $F1"
    $repeatscout_path/RepeatScout -sequence $assembly_dir/$F1/$F1.fasta -output $assembly_dir/$F1/repeatscout/$F1.output_repeats.fas -freq $assembly_dir/$F1/repeatscout/$F1.output_lmer.frequency
    echo "running step-3:filter-stage-1.prl for $F1"
    $repeatscout_path/filter-stage-1.prl $assembly_dir/$F1/repeatscout/$F1.output_repeats.fas > $assembly_dir/$F1/repeatscout/$F1.output_repeats.fas.filtered_1 
    echo "running step-4:filter-stage-2.prl for $F1"
    $repeatmasker_path/RepeatMasker $assembly_dir/$F1/$F1.fasta -lib $assembly_dir/$F1/repeatscout/$F1.output_repeats.fas.filtered_1 -dir $assembly_dir/$F1/repeatmasker/ -no_is
    echo "running step-5:filter-stage-2.prl for $F1"
    cat $assembly_dir/$F1/repeatscout/$F1.output_repeats.fas.filtered_1 | $repeatscout_path/filter-stage-2.prl --cat=$assembly_dir/$F1/repeatmasker/$F1.fasta.out > $assembly_dir/$F1/repeatscout/output_repeats.fas.filtered_2
    #cat repeats.lib | ./filter-stage-2.prl --cat=repeats.out --thresh=20

    cp $assembly_dir/$F1/repeatmasker/$F1.fasta.masked $assembly_dir/all_fasta/$F1.repeatscout_masked.fasta
done

###############################################################################
echo "completed.... repeatscout and RepeatMasker ------------------------------"
###############################################################################