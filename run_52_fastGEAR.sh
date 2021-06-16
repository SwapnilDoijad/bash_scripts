#!/bin/bash
###############################################################################
## Note: provide full path, otherwise sometime doesenot work
###############################################################################
echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

echo "want to suffix files? type "y" else press enter to continue"
read answer
if [ "$answer" = "y" ]; then
s="_$l"
fi

###############################################################################

mkdir results/52_fastGEAR_ancestral"$s"

bash /home/swapnil/tools/fastGEARpackageLinux64bit/run_fastGEAR.sh /usr/local/MATLAB/MATLAB_Runtime/v901 results/17_roary"$s"/roary_results/core_gene_alignment.aln results/52_fastGEAR_ancestral"$s"/out.mat /home/swapnil/tools/fastGEARpackageLinux64bit/fG_input_specs.txt

###############################################################################
