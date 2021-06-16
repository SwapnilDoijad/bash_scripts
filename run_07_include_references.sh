#!/bin/bash
###############################################################################
# including references
###############################################################################

echo "started... step-7 including references ---------------------------------"

###############################################################################

if [ -d data/references/ ]; then
cp data/references/*.fasta results/00_ref/
fi

cp results/00_ref/*.fasta results/04_assembly/all_fasta/

if [ -f results/04_assembly/all_fasta/ref.*.fasta ]; then
rename 's/ref\.//' results/04_assembly/all_fasta/ref.*.fasta
fi

cd results/00_ref/
ls *.fasta > tmp.tmp
sed -i 's/\.fasta//g' tmp.tmp
cd ..
cd ..

cp list.txt list_wo_ref.txt

cat results/00_ref/tmp.tmp >> list.txt
sed -i 's/ref\.//g' list.txt
sed -i '/^\s*$/d' list.txt

rm results/00_ref/tmp.tmp


###############################################################################

echo "completed... step-7 including references -------------------------------"

###############################################################################
