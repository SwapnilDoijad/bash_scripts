#!/bin/bash
###############################################################################
#31 WGS-alignment-mugsy

###############################################################################

echo "started.... step-31 mugsy ----------------------------------------------"

###############################################################################
if [ ! -d results/31_mugsy ] ; then
mkdir results/31_mugsy
fi

if [ ! -d results/31_mugsy/fasta ] ; then
mkdir results/31_mugsy/fasta
fi

#------------------------------------------------------------------------------

#cp results/0_ref/ref.*.fasta results/31_mugsy/fasta/
#cp results/2_assembly/all_fasta/*.fasta results/31_mugsy/fasta/

#------------------------------------------------------------------------------
exit

source /home/swapnil/tools/mugsy/mugsyenv.sh
(/home/swapnil/pipeline/tools/mugsy/mugsy --directory /home/swapnil/pipeline/results/31_mugsy --prefix mugsy_alignment /home/swapnil/pipeline/results/31_mugsy/fasta/*.fasta) > /dev/null 2>&1

V1=$(ls -1 results/31_mugsy/fasta | wc -l)
awk -f /home/swapnil/pipeline/tools/maf2phy.awk -v n=$V1 results/31_mugsy/mugsy_alignment.maf > results/31_mugsy/mugsy_alignment.maf.phy

cp results/31_mugsy/mugsy_alignment.maf.phy /home/swapnil/pipeline/tools/
cd /home/swapnil/pipeline/tools
python2.7 convert_sequence_format_for_run_31_mugsy.py
cd ..
rm /home/swapnil/pipeline/tools/mugsy_alignment.maf.phy
cp /home/swapnil/pipeline/tools/mugsy_alignment.maf.phy.fa results/31_mugsy/


(raxmlHPC -m GTRGAMMA -p 12345 -x 12345 -# 100 -s results/31_mugsy/mugsy_alignment.maf.phy.fa -n T19) > /dev/null 2>&1
(raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s results/31_mugsy/mugsy_alignment.maf.phy.fa -n T20) > /dev/null 2>&1

mv RAxML_bipartitions.T20 results/31_mugsy/RAxML_bipartitions.T20.NO-Rec-Filt.tree
rm *.T19
rm *.T20

#------------------------------------------------------------------------------


if [ ! -d results/31_mugsy/gubbins ] ; then
mkdir results/31_mugsy/gubbins
fi

(run_gubbins results/31_mugsy/mugsy_alignment.maf.phy.fa) > /dev/null 2>&1

mv mugsy_alignment.* results/31_mugsy/gubbins
mv *.mugsy_alignment.* results/31_mugsy/gubbins

(raxmlHPC -m GTRGAMMA -p 12345 -x 12345 -# 100 -s results/31_mugsy/gubbins/mugsy_alignment.maf.phy.filtered_polymorphic_sites.fasta -n T19) > /dev/null 2>&1
(raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s results/31_mugsy/gubbins/mugsy_alignment.maf.phy.filtered_polymorphic_sites.fasta -n T20) > /dev/null 2>&1
mv RAxML_bipartitions.T20 results/31_mugsy/gubbins/RAxML_bipartitions.T20.Rec-region-removed.tree
rm *.T19
rm *.T20

###############################################################################

echo "completed.... step-31 mugsy --------------------------------------------"

###############################################################################
