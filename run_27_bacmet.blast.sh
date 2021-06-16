#27 BacMet

###############################################################################
echo "started.... BacMet -----------------------------------------------------"
###############################################################################

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | sed 's/list\.//g' | sed 's/\.txt//g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")
	s="_$l"

	(mkdir results/27_BacMet"$s") > /dev/null 2>&1
	(mkdir results/27_BacMet"$s"/tmp) > /dev/null 2>&1
	(mkdir results/27_BacMet"$s"/tmp_EXP) > /dev/null 2>&1
	(mkdir results/27_BacMet"$s"/tmp_PRE) > /dev/null 2>&1

###############################################################################
for F1 in $(cat $list); do
	## chekc if already finished
	if [ -f results/27_BacMet"$s"/tmp_EXP/$F1.BacMet_EXP.tmp ]; then
		echo "bacmet analysis for $F1 already finished"
		else

		echo "runnning... BacMet for $F1"

		blastp -db /media/swapnil/share/databases/bacmet/BacMet_EXP.704.Header_only_Bac-Id.fasta -query results/08_annotation/raw_files/$F1/$F1.faa -out results/27_BacMet"$s"/tmp_EXP/$F1.BacMet_EXP.tmp -max_target_seqs 1 -evalue 1e-100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"

		blastp -db /media/swapnil/share/databases/bacmet/BacMet_PRE.40556.Header_only_Gi.fasta -query results/08_annotation/raw_files/$F1/$F1.faa -out results/27_BacMet"$s"/tmp_PRE/$F1.BacMet_PRE.tmp -max_target_seqs 1 -evalue 1e-100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"

		awk '{print $2}' results/27_BacMet"$s"/tmp_EXP/$F1.BacMet_EXP.tmp > results/27_BacMet"$s"/tmp/$F1.BacMet_EXP.list.tmp
			
			for F2 in $(cat results/27_BacMet"$s"/tmp/$F1.BacMet_EXP.list.tmp);do
			sed -i '/'$F2'/d' results/27_BacMet"$s"/tmp_PRE/$F1.BacMet_PRE.tmp
			done

		cat results/27_BacMet"$s"/tmp_EXP/$F1.BacMet_EXP.tmp results/27_BacMet"$s"/tmp_PRE/$F1.BacMet_PRE.tmp > results/27_BacMet"$s"/tmp/$F1.BacMet.tmp

		#------------------------------------------------------------------------------
		awk '{print $1}' results/27_BacMet"$s"/tmp/$F1.BacMet.tmp > results/27_BacMet"$s"/tmp/$F1.BacMet.list.tmp

		awk -F'\t' -f /home/swapnil/pipeline/tools/vlookup-BacMet.awk /media/swapnil/share/databases/bacmet/BacMet_both.mapping.tab results/27_BacMet"$s"/tmp/$F1.BacMet.list.tmp > results/27_BacMet"$s"/tmp/$F1.BacMet.list.metadata.tmp

		paste results/27_BacMet"$s"/tmp/$F1.BacMet.tmp results/27_BacMet"$s"/tmp/$F1.BacMet.list.metadata.tmp > results/27_BacMet"$s"/$F1.BacMet.tab

		cat results/27_BacMet"$s"/$F1.BacMet.tab >> results/27_BacMet"$s"/all.BacMet.tab

		ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs gene protein group originated-from gene-name location organism compound description' -cx results/27_BacMet"$s"/$F1.BacMet.tab
		sed -i '1s/ /\t/g' results/27_BacMet"$s"/$F1.BacMet.tab

		#------------------------------------------------------------------------------

		echo "runnning... BacMet for $F1"

	fi
done

ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs gene protein group originated-from gene-name location organism compound description' -cx results/27_BacMet"$s"/all.BacMet.tab
sed -i '1s/ /\t/g' results/27_BacMet"$s"/all.BacMet.tab
###############################################################################
## Creating matrix

echo "creating matrix for BacMet"

for F1 in $(cat $list);do

	awk 'FNR>1 {print $16}' results/27_BacMet"$s"/all.BacMet.tab | sort | uniq > results/27_BacMet"$s"/tmp/all.BacMet.list.tab

	awk 'FNR>1 {print $16}' results/27_BacMet"$s"/$F1.BacMet.tab > results/27_BacMet"$s"/tmp/$F1.BacMet.gene-name-list.tab

		for F2 in $(cat results/27_BacMet"$s"/tmp/all.BacMet.list.tab);do
		V2=$(awk '{count[$1]++} END {print count["'$F2'"]}' results/27_BacMet"$s"/tmp/$F1.BacMet.gene-name-list.tab)
		echo $V2 >> results/27_BacMet"$s"/tmp/$F1.BacMet.list.count.1.tab
		done
	awk '{for (i=1; i<= NF; i++) {if($i > 1) { $i=1; } } print }' results/27_BacMet"$s"/tmp/$F1.BacMet.list.count.1.tab > results/27_BacMet"$s"/tmp/$F1.BacMet.list.count.2.tab
	ex -sc '1i|'$F1'' -cx results/27_BacMet"$s"/tmp/$F1.BacMet.list.count.1.tab
	ex -sc '1i|'$F1'' -cx results/27_BacMet"$s"/tmp/$F1.BacMet.list.count.2.tab
	sed -i -e 's/^$/0/' results/27_BacMet"$s"/tmp/$F1.BacMet.list.count.1.tab
	sed -i -e 's/^$/0/' results/27_BacMet"$s"/tmp/$F1.BacMet.list.count.2.tab

	cp results/27_BacMet"$s"/tmp/all.BacMet.list.tab results/27_BacMet"$s"/tmp/all.BacMet.list.2.tab
	sed -i 's/$/:c/' results/27_BacMet"$s"/tmp/all.BacMet.list.2.tab
	sed -i '1i\===\' results/27_BacMet"$s"/tmp/all.BacMet.list.2.tab
	paste results/27_BacMet"$s"/tmp/all.BacMet.list.2.tab results/27_BacMet"$s"/tmp/*.BacMet.list.count.1.tab > results/27_BacMet"$s"/tmp/matrix.csv.tmp
	paste results/27_BacMet"$s"/tmp/all.BacMet.list.2.tab results/27_BacMet"$s"/tmp/*.BacMet.list.count.2.tab > results/27_BacMet"$s"/tmp/matrix_1-0.csv.tmp

		##-----------------------------------------------------------------------------
		## transpose
		awk '
		{ 
			for (i=1; i<=NF; i++)  {
				a[NR,i] = $i
			}
		}
		NF>p { p = NF }
		END {    
			for(j=1; j<=p; j++) {
				str=a[1,j]
				for(i=2; i<=NR; i++){
					str=str" "a[i,j];
				}
				print str
			}
		}' results/27_BacMet"$s"/tmp/matrix.csv.tmp > results/27_BacMet"$s"/matrix.csv

		awk '
		{ 
			for (i=1; i<=NF; i++)  {
				a[NR,i] = $i
			}
		}
		NF>p { p = NF }
		END {    
			for(j=1; j<=p; j++) {
				str=a[1,j]
				for(i=2; i<=NR; i++){
					str=str" "a[i,j];
				}
				print str
			}
		}' results/27_BacMet"$s"/tmp/matrix_1-0.csv.tmp > results/27_BacMet"$s"/matrix_1-0.csv
		##-----------------------------------------------------------------------------

	sed -i 's/===//g' results/27_BacMet"$s"/matrix.csv
	sed -i 's/===//g' results/27_BacMet"$s"/matrix_1-0.csv

done

echo "Completed matrix for BacMet"
###############################################################################
## Count frequency
for F1 in $(cat $list);do
awk 'FNR>1 {print $16}' results/27_BacMet"$s"/$F1.BacMet.tab >> results/27_BacMet"$s"/tmp/all.BacMet.gene-list.tab
done
	for F2 in $(cat results/27_BacMet"$s"/tmp/all.BacMet.list.tab);do
	V2=$(awk '{count[$1]++} END {print count["'$F2'"]}' results/27_BacMet"$s"/tmp/all.BacMet.gene-list.tab)
	echo $F2 $V2 >> results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.tmp
	done

cat results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.tmp | sort -rnk 2 -t' ' > results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.sorted.tmp
sed -i 's/ /\t/g' results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.sorted.tmp
rm results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.tmp



awk -F'\t' -f /home/swapnil/pipeline/tools/vlookup-BacMet-YIN_new.awk /media/swapnil/share/databases/bacmet/BacMet_both.mapping.2.tab results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.sorted.tmp > results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.sorted.metadata.tmp

paste results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.sorted.tmp results/27_BacMet"$s"/tmp/all.BacMet.gene-list.count.sorted.metadata.tmp > results/27_BacMet"$s"/all.BacMet.count.tab

ex -sc '1i|gene frequency location organism compound description' -cx results/27_BacMet"$s"/all.BacMet.count.tab
sed -i '1s/ /\t/g' results/27_BacMet"$s"/all.BacMet.count.tab

unoconv -i FilterOptions=09,,system,1 -f xls results/27_BacMet"$s"/all.BacMet.count.tab


###############################################################################

echo "Completed.... BacMet ---------------------------------------------------"

###############################################################################
