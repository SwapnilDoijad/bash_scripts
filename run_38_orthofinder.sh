#!/bin/bash
###############################################################################
#XX orthofinder
###############################################################################

echo "started.... 38_OrthoFinder ---------------------------------------------"

###############################################################################
if [ ! -d results/38_OrthoFinder ]; then
mkdir results/38_OrthoFinder
fi

tools/OrthoFinder-2.2.3/orthofinder -f results/08_annotation/all_faa/


###############################################################################

echo "completed.... 38_OrthoFinder -------------------------------------------"

###########################################################################
