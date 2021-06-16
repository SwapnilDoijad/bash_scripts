source /home/swapnil/miniconda3/etc/profile.d/conda.sh
conda activate myenv
###############################################################################
for F1 in $(cat list.txt);do
mkdir $F1
# falco /media/swapnil/network/$F1.fastq.gz  -d /tmp -o $F1/
fastqc /media/swapnil/network/$F1.fastq.gz
mulitqc 
wkhtmltopdf $F1/fastqc_report.html $F1/$F1.pdf
done
###############################################################################
conda deactivate