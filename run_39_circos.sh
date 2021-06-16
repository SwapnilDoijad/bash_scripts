#!/bin/bash
###############################################################################
#35 Circos
###############################################################################
#echo "Creating Circos images--------------------------------------------------"
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

	echo "data to incldue?"
	echo "to include virulence gene? type y"
	read virulence_genes
		if [ "$virulence_genes" = "y" ]; then
		echo "provide virulence gene folder path (for e.g. results/21_VFDB_blast_results)"
		read VG_path
		fi

	echo "to include antibiotic resistance gene? type y"
	read Abr_genes
		if [ "$Abr_genes" = "y" ]; then
		echo "provide Abr gene folder path (for e.g. results/22_CARD-AR_results/results)"
		read abr_path
		fi

	echo "OK"
	(mkdir results/39_circos) > /dev/null 2>&1
	(mkdir results/39_circos/results) > /dev/null 2>&1
	(mkdir results/39_circos/tmp) > /dev/null 2>&1

	path=$(pwd)
###############################################################################
#collect data

for F1 in $(cat $list); do
	if [ ! -f results/39_circos/results/$F1/$F1.circos.png ]; then 

		(mkdir results/39_circos/results/$F1) > /dev/null 2>&1
		(mkdir results/39_circos/tmp/$F1) > /dev/null 2>&1
		#genome size
		V1=$(awk 'FNR == 1 {print $2}' results/08_annotation/raw_files/$F1/tmp/$F1.08_annotation.statistics.txt)

		#tRNA
		cat results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' | awk -F',' '$1 == "tRNA" {print ($3-4000) " " ($4+4000)}'  > results/39_circos/tmp/$F1/tRNA.txt
		sed -i 's/-4000//' results/39_circos/tmp/$F1/tRNA.txt
		sed -i 's/ 4000//' results/39_circos/tmp/$F1/tRNA.txt
		sed -i 's/^/chr1 /' results/39_circos/tmp/$F1/tRNA.txt
		sed -i 's/$/ color=255,20,147/' results/39_circos/tmp/$F1/tRNA.txt

		#rRNA
		cat results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' | awk -F',' '$1 == "rRNA" {print ($3-3000) " " ($4+3000)}'  > results/39_circos/tmp/$F1/rRNA.txt
		sed -i 's/-3000//' results/39_circos/tmp/$F1/tRNA.txt
		sed -i 's/ 3000//' results/39_circos/tmp/$F1/tRNA.txt
		sed -i 's/^/chr1 /' results/39_circos/tmp/$F1/rRNA.txt
		sed -i 's/$/ color=0,100,0/' results/39_circos/tmp/$F1/rRNA.txt

		#CDS-F
		cat results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' | awk -F',' '($1 == "CDS" && $6 == "no") {print $3,$4}' > results/39_circos/tmp/$F1/CDS_F.txt
		sed -i 's/^/chr1 /' results/39_circos/tmp/$F1/CDS_F.txt
		sed -i 's/$/ color=0,255,0/' results/39_circos/tmp/$F1/CDS_F.txt

		#CDS-R
		cat results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' | awk -F',' '($1 == "CDS" && $6 == "yes") {print $3,$4}' > results/39_circos/tmp/$F1/CDS_R.txt
		sed -i 's/^/chr1 /' results/39_circos/tmp/$F1/CDS_R.txt
		sed -i 's/$/ color=255,140,0/' results/39_circos/tmp/$F1/CDS_R.txt

		#Virulence-genes					
		if [ "$virulence_genes" = "y" ]; then
			VG_file=$VG_path/alignment_results/$F1.VFDB-blast.csv
			awk 'NR>1 {print $2}' $VG_file > results/39_circos/tmp/$F1/VG.tmp.1.txt
				(rm results/39_circos/tmp/$F1/VG.tmp.2.txt) > /dev/null 2>&1
				locus_tag=$(head -1 results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' | awk -F',' '{ for (i=1;i<=NF;i++) if ($i == "locus_tag") print i }')
				for F2 in $(cat results/39_circos/tmp/$F1/VG.tmp.1.txt); do
				vg_coordinates=$(cat results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' | awk -F',' '( $1 == "CDS" && $"'$locus_tag'" == "'$F2'") {print $3,$4}')
                echo $vg_coordinates >> results/39_circos/tmp/$F1/VG.tmp.2.txt
                done
			sed 's/^/chr1 /' results/39_circos/tmp/$F1/VG.tmp.2.txt > results/39_circos/tmp/$F1/VG.txt
			sed -i 's/$/ color=dred/' results/39_circos/tmp/$F1/VG.txt

			#Virulence-genes-label
			(rm results/39_circos/tmp/$F1/VG.label.txt) > /dev/null 2>&1
			gene=$(head -1 $VG_file | sed 's/"//g' | awk -F'\t' '{ for (i=1;i<=NF;i++) if ($i == "gene") print i }')
			awk -F'\t' 'NR>1 {print $'$gene' }' $VG_file | sed 's/ /_/g' > results/39_circos/tmp/$F1/VG.tmp.3.txt
			paste results/39_circos/tmp/$F1/VG.tmp.2.txt results/39_circos/tmp/$F1/VG.tmp.3.txt > results/39_circos/tmp/$F1/VG.label.txt
			sed -i 's/^/chr1 /' results/39_circos/tmp/$F1/VG.label.txt
			sed -i 's/\t/ /g' results/39_circos/tmp/$F1/VG.label.txt

			# count number of VG labels an fix maximum to 20
			VG_total=$(cat results/39_circos/tmp/$F1/VG.label.txt | wc -l)

			if [ $VG_total -gt 15 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/VG.label.txt ; fi
			VG_total=$(cat results/39_circos/tmp/$F1/VG.label.txt | wc -l)

			if [ $VG_total -gt 15 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/VG.label.txt ; fi			
			VG_total=$(cat results/39_circos/tmp/$F1/VG.label.txt | wc -l)

			if [ $VG_total -gt 15 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/VG.label.txt ; fi
			VG_total=$(cat results/39_circos/tmp/$F1/VG.label.txt | wc -l)

			if [ $VG_total -gt 15 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/VG.label.txt ; fi			
			VG_total=$(cat results/39_circos/tmp/$F1/VG.label.txt | wc -l)

			if [ $VG_total -gt 15 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/VG.label.txt ; fi
		fi

		## Abr-genes
		if [ "$Abr_genes" = "y" ]; then
			abr_file=$abr_path/$F1.txt
			awk -F' ' 'NR>1 {print $1}' $abr_file > results/39_circos/tmp/$F1/Abr.tmp.1.txt
				(rm results/39_circos/tmp/$F1/Abr.tmp.2.txt) > /dev/null 2>&1
				locus_tag=$(head -1 results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' | awk -F',' '{ for (i=1;i<=NF;i++) if ($i == "locus_tag") print i }')
				for F2 in $(cat results/39_circos/tmp/$F1/Abr.tmp.1.txt); do
				Abr_coordinates=$(cat results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' | awk -F',' '( $1 == "CDS" && $"'$locus_tag'" == "'$F2'") {print $3,$4}')
				echo $Abr_coordinates >> results/39_circos/tmp/$F1/Abr.tmp.2.txt
				done
			sed 's/^/chr1 /' results/39_circos/tmp/$F1/Abr.tmp.2.txt > results/39_circos/tmp/$F1/Abr.txt
			sed -i 's/$/ color=dred/' results/39_circos/tmp/$F1/Abr.txt

			#Abr-genes-label
			(rm results/39_circos/tmp/$F1/Abr.tmp.3a.txt) > /dev/null 2>&1

			awk -F'\t' 'NR>1 {print $9}' $abr_file | sed "s/ /_/g"  > results/39_circos/tmp/$F1/Abr.tmp.3.txt
				(rm results/39_circos/tmp/$F1/Abr.tmp.3a.txt) > /dev/null 2>&1
				for F3 in $(cat results/39_circos/tmp/$F1/Abr.tmp.3.txt); do
				lable_words=$(echo $F3 | wc -c )
				if [ "$lable_words" -gt 15 ] ; then
				(echo $F3 | awk -F'_' '{print $3}') >> results/39_circos/tmp/$F1/Abr.tmp.3a.txt
				else
				echo $F3 >> results/39_circos/tmp/$F1/Abr.tmp.3a.txt
				fi
				done
			paste results/39_circos/tmp/$F1/Abr.tmp.2.txt results/39_circos/tmp/$F1/Abr.tmp.3a.txt > results/39_circos/tmp/$F1/Abr.label.txt
			sed -i 's/^/chr1 /' results/39_circos/tmp/$F1/Abr.label.txt
			sed -i 's/\t/ /g' results/39_circos/tmp/$F1/Abr.label.txt

			# count number of VG labels an fix maximum to 20
			ABR_total=$(cat results/39_circos/tmp/$F1/Abr.label.txt| wc -l)			
			if [ $ABR_total -gt 20 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/Abr.label.txt; fi

			ABR_total=$(cat results/39_circos/tmp/$F1/Abr.label.txt| wc -l)			
			if [ $ABR_total -gt 20 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/Abr.label.txt; fi

			ABR_total=$(cat results/39_circos/tmp/$F1/Abr.label.txt| wc -l)			
			if [ $ABR_total -gt 20 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/Abr.label.txt; fi

			ABR_total=$(cat results/39_circos/tmp/$F1/Abr.label.txt| wc -l)			
			if [ $ABR_total -gt 20 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/Abr.label.txt; fi

			ABR_total=$(cat results/39_circos/tmp/$F1/Abr.label.txt| wc -l)			
			if [ $ABR_total -gt 20 ] ; then  sed -i '1~3d' results/39_circos/tmp/$F1/Abr.label.txt; fi
		fi

		## custom gene set-1
		if [ "$custom_gene_set_1" = "y" ]; then
			custom_gene_set_1_file=$(find $custom_gene_set_1_path/$F1* )
			echo $custom_gene_set_1_file
			#custom gene set-1 co-ordinates and color
			awk 'NR==1{for(i=1; i<=NF; i++) if ( $i == "qstart" ) {a[i]++;} } { for (i in a) printf "%s\t", ($i-10000); printf "\n"}' $custom_gene_set_1_file  | sed '1d' > results/39_circos/tmp/$F1/custom_gene_set_1.start.tmp
			awk 'NR==1{for(i=1; i<=NF; i++) if ( $i == "qend" ) {a[i]++;} } { for (i in a) printf "%s\t", ($i+10000); printf "\n"}' $custom_gene_set_1_file | sed '1d' > results/39_circos/tmp/$F1/custom_gene_set_1.end.tmp
			paste results/39_circos/tmp/$F1/custom_gene_set_1.start.tmp results/39_circos/tmp/$F1/custom_gene_set_1.end.tmp > results/39_circos/tmp/$F1/custom_gene_set_1.coordinates.txt
			sed -i 's/$/ color=dred/' results/39_circos/tmp/$F1/custom_gene_set_1.coordinates.txt
			sed -i 's/^/chr1 /' results/39_circos/tmp/$F1/custom_gene_set_1.coordinates.txt
			#custom gene set-1 label
			awk 'NR==1{for(i=1; i<=NF; i++) if ($i=="qstart" || $i == "qend" ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' $custom_gene_set_1_file | sed 's/^/chr1 /' | sed '1d' > results/39_circos/tmp/$F1/custom_gene_set_1.labels.1.tmp
			awk 'NR==1{for(i=1; i<=NF; i++) if ($i == "sseqid" ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' $custom_gene_set_1_file | sed '1d'  > results/39_circos/tmp/$F1/custom_gene_set_1.labels.2.tmp
			paste results/39_circos/tmp/$F1/custom_gene_set_1.labels.1.tmp results/39_circos/tmp/$F1/custom_gene_set_1.labels.2.tmp | sed 's/\t/ /g' > results/39_circos/tmp/$F1/custom_gene_set_1.labels.txt
			rm results/39_circos/tmp/$F1/*.tmp
		fi

		###############################################################################
		#write collected all data

		echo "chr - chr1 1 0 $V1 black" > results/39_circos/tmp/$F1/genome_size.txt
		echo "chr1 0 $V1" > results/39_circos/tmp/$F1/highlight.txt
		echo "chr1 0 $V1" > results/39_circos/tmp/$F1/seperator.txt

		#------------------------------------------------------------------------------
		#Copy all data

		cp results/39_circos/tmp/$F1/*.txt /home/swapnil/tools/circos-0.69/
		cp results/08_annotation/raw_files/$F1/$F1.fna /home/swapnil/tools/circos-0.69/$F1.genome.txt

		#------------------------------------------------------------------------------
		#run circos

		cd /home/swapnil/tools/circos-0.69

		(perl GCcalc-master/gcSkew.pl -f $F1.genome.txt -o $F1.gcskew.txt) > /dev/null 2>&1
		awk -F'\t' '{print $1,$2,$3,$4}' $F1.gcskew.txt > gcskew.txt

		echo "running circos for $F1"
		#(
		./bin/circos -conf conf/conf.VG.Abr.conf
		#) > /dev/null 2>&1

		cd $path

		#------------------------------------------------------------------------------
		#move image and circos-files from home directory to result directory

		mv /home/swapnil/tools/circos-0.69/circos.png results/39_circos/results/$F1/$F1.circos.png ;
		mv /home/swapnil/tools/circos-0.69/circos.svg results/39_circos/results/$F1/$F1.circos.svg ;
		rm /home/swapnil/tools/circos-0.69/*.txt

		#------------------------------------------------------------------------------
		else
		echo "$F1 is already finished"
	fi
done
###############################################################################
#echo "Finisehd creating Circos images-----------------------------------------"
###############################################################################

