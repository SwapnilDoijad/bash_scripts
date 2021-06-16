#!/bin/bash
###############################################################################
#28 traitar

###############################################################################

echo "started.... traitar ----------------------------------------------------"

###############################################################################
if [ ! -d results/28_traitar ]; then
mkdir results/28_traitar
fi

if [ -f data/traitar_samples.txt ]; then
	cp data/traitar_samples.txt results/28_traitar/samples.txt
else
	cp list.txt results/28_traitar/tmp1.tmp
	sed -i 's/$/\.faa/' results/28_traitar/tmp1.tmp
	cp list.txt results/28_traitar/tmp2.tmp
	sed -i -e 's/^/strain_/' results/28_traitar/tmp2.tmp
	sed -i 's/$/	none/' results/28_traitar/tmp2.tmp

	paste results/28_traitar/tmp1.tmp results/28_traitar/tmp2.tmp > results/28_traitar/samples.txt

	ex -sc '1i|sample_file_name sample_name category' -cx results/28_traitar/samples.txt
	sed -i '1s/ /\t/g' results/28_traitar/samples.txt

	rm results/28_traitar/*.tmp
fi

traitar phenotype results/08_annotation/all_faa results/28_traitar/samples.txt from_genes results/28_traitar/results

###############################################################################

echo "Completed.... traitar --------------------------------------------------"

###############################################################################




