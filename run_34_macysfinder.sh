#!/bin/bash
###############################################################################
# 34 MacSyFinder
###############################################################################
echo "started... step-34 MacSyFinder -----------------------------------------"
###############################################################################

echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")


(mkdir results/34_MacSyFinder) > /dev/null 2>&1
(mkdir results/34_MacSyFinder/results_strict) > /dev/null 2>&1
(mkdir results/34_MacSyFinder/results_strict/report) > /dev/null 2>&1
(mkdir results/34_MacSyFinder/results_loose) > /dev/null 2>&1
(mkdir results/34_MacSyFinder/results_loose/report) > /dev/null 2>&1

for F1 in $(cat $list); do

echo "running..... step-34 MacSyFinder for..." $F1

#-----------------------
## system strict cutoff
(mkdir results/34_MacSyFinder/results_strict/$F1) > /dev/null 2>&1
(mkdir results/34_MacSyFinder/results_strict/$F1/tmp) > /dev/null 2>&1

echo "#Replicon_name	System_Id	Reference_system	System_status	Nb_loci	Nb_Ref_mandatory	Nb_Ref_accessory	Nb_Ref_Genes_detected_NR	Nb_Genes_with_match	System_length	Nb_Mandatory_NR	Nb_Accessory_NR	Nb_missing_mandatory	Nb_missing_accessory	List_missing_mandatory	List_missing_accessory	Loci_positions	Occur_Mandatory	Occur_Accessory	Occur_Forbidden" > results/34_MacSyFinder/results_strict/report/$F1.summary.tab

echo "#Hit_Id	Replicon_name	Position	Sequence_length	Gene	Reference_system	Predicted_system	System_Id	System_status	Gene_status	i-evalue	Score	Profile_coverage	Sequence_coverage	Begin_match	End_match" > results/34_MacSyFinder/results_strict/report/$F1.report.tab

	for F2 in $(cat /media/swapnil/share/databases/macsyfinder/system.list.txt);do

	(macsyfinder -d /media/swapnil/share/databases/macsyfinder/DEF/  $F2 -p /media/swapnil/share/databases/macsyfinder/profiles/ -o results/34_MacSyFinder/results_strict/$F1/$F1.$F2 --db-type ordered_replicon --sequence-db results/08_annotation/raw_files/$F1/$F1.faa) > /dev/null 2>&1
	
	rm -rf results/34_MacSyFinder/results_strict/$F1/$F1.$F2/hmmer_results
	rm results/34_MacSyFinder/results_strict/$F1/$F1.$F2/macsyfinder.conf
	rm results/34_MacSyFinder/results_strict/$F1/$F1.$F2/macsyfinder.log
	rm results/34_MacSyFinder/results_strict/$F1/$F1.$F2/macsyfinder.out

	if [[ -s results/34_MacSyFinder/results_strict/$F1/$F1.$F2/macsyfinder.summary ]] ; then
	awk 'FNR == 2 {print}' results/34_MacSyFinder/results_strict/$F1/$F1.$F2/macsyfinder.summary >> results/34_MacSyFinder/results_strict/report/$F1.summary.tab
	awk 'FNR >1 {print}' results/34_MacSyFinder/results_strict/$F1/$F1.$F2/macsyfinder.report >> results/34_MacSyFinder/results_strict/report/$F1.report.tab
	else
	:
	fi ;

	done
#-----------------------
awk -F'\t' 'FNR >1 {print $3}' results/34_MacSyFinder/results_strict/report/$F1.summary.tab > results/34_MacSyFinder/results_strict/$F1/tmp/$F1.summary.tab.tmp
echo $F1 > results/34_MacSyFinder/results_strict/$F1/tmp/$F1.summary.tab.2.tmp

for F2 in $(cat /media/swapnil/share/databases/macsyfinder/system.list.txt);do
grep -c "$F2" results/34_MacSyFinder/results_strict/$F1/tmp/$F1.summary.tab.tmp >> results/34_MacSyFinder/results_strict/$F1/tmp/$F1.summary.tab.2.tmp ;
done

#-----------------------
## system loose cutoff

mkdir results/34_MacSyFinder/results_loose/$F1
mkdir results/34_MacSyFinder/results_loose/$F1/tmp

echo "#Replicon_name	System_Id	Reference_system	System_status	Nb_loci	Nb_Ref_mandatory	Nb_Ref_accessory	Nb_Ref_Genes_detected_NR	Nb_Genes_with_match	System_length	Nb_Mandatory_NR	Nb_Accessory_NR	Nb_missing_mandatory	Nb_missing_accessory	List_missing_mandatory	List_missing_accessory	Loci_positions	Occur_Mandatory	Occur_Accessory	Occur_Forbidden" > results/34_MacSyFinder/results_loose/report/$F1.summary.tab

echo "#Hit_Id	Replicon_name	Position	Sequence_length	Gene	Reference_system	Predicted_system	System_Id	System_status	Gene_status	i-evalue	Score	Profile_coverage	Sequence_coverage	Begin_match	End_match" > results/34_MacSyFinder/results_loose/report/$F1.report.tab

	for F2 in $(cat /media/swapnil/share/databases/macsyfinder/system.list.txt);do

	(macsyfinder -d /media/swapnil/share/databases/macsyfinder/DEF/  $F2 -p /media/swapnil/share/databases/macsyfinder/profiles/ -o results/34_MacSyFinder/results_loose/$F1/$F1.$F2 --db-type unordered_replicon --sequence-db results/08_annotation/raw_files/$F1/$F1.faa) > /dev/null 2>&1
	

	rm -rf results/34_MacSyFinder/results_loose/$F1/$F1.$F2/hmmer_results
	rm results/34_MacSyFinder/results_loose/$F1/$F1.$F2/macsyfinder.conf
	rm results/34_MacSyFinder/results_loose/$F1/$F1.$F2/macsyfinder.log
	rm results/34_MacSyFinder/results_loose/$F1/$F1.$F2/macsyfinder.out

	if [[ -s results/34_MacSyFinder/results_loose/$F1/$F1.$F2/macsyfinder.summary ]] ; then
	awk 'FNR == 2 {print}' results/34_MacSyFinder/results_loose/$F1/$F1.$F2/macsyfinder.summary >> results/34_MacSyFinder/results_loose/report/$F1.summary.tab
	awk 'FNR >1 {print}' results/34_MacSyFinder/results_loose/$F1/$F1.$F2/macsyfinder.report >> results/34_MacSyFinder/results_loose/report/$F1.report.tab
	else
	:
	fi ;

	done
#---------------------	

awk -F'\t' 'FNR >1 {print $3}' results/34_MacSyFinder/results_loose/report/$F1.summary.tab > results/34_MacSyFinder/results_loose/$F1/tmp/$F1.summary.tab.tmp
echo $F1 > results/34_MacSyFinder/results_loose/$F1/tmp/$F1.summary.tab.2.tmp

for F2 in $(cat /media/swapnil/share/databases/macsyfinder/system.list.txt);do
grep -c "$F2" results/34_MacSyFinder/results_loose/$F1/tmp/$F1.summary.tab.tmp >> results/34_MacSyFinder/results_loose/$F1/tmp/$F1.summary.tab.2.tmp ;
done

#----------
echo "finished.... step-34 MacSyFinder for..." $F1

done

paste /media/swapnil/share/databases/macsyfinder/system.list.2.txt results/34_MacSyFinder/results_strict/*/tmp/*.summary.tab.2.tmp > results/34_MacSyFinder/results_strict/all.tab.tmp
paste /media/swapnil/share/databases/macsyfinder/system.list.2.txt results/34_MacSyFinder/results_loose/*/tmp/*.summary.tab.2.tmp > results/34_MacSyFinder/results_loose/all.tab.tmp

sh /home/swapnil/pipeline/tools/transpose-macsyfinder.sh

rm results/34_MacSyFinder/results_strict/all.tab.tmp
rm results/34_MacSyFinder/results_loose/all.tab.tmp

###############################################################################
echo "completed... step-34 MacSyFinder ---------------------------------------"
###############################################################################
