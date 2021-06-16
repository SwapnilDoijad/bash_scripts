#!/bin/bash
###############################################################################
#25 plasmid serach complete
## June 2019 needs to refine ? 
###############################################################################
echo "started.... plasmid search complete -------------------------------------"
###############################################################################
## preliminary file preparation
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")
    s=$(echo "_$l")

	#------------------------------------------------------------------------------
	(mkdir results/25_plasmid_finder_complete_BLAST"$s") > /dev/null 2>&1
	(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files) > /dev/null 2>&1
###############################################################################
## search plasmid

	for F1 in $(cat $list); do
			echo "runnning... Plasmid search for $F1 .........................."
			(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1) > /dev/null 2>&1
			(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS) > /dev/null 2>&1
			(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp) > /dev/null 2>&1
			(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta) > /dev/null 2>&1
			(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta_RefSeq_plasmid_first_attmpt_shortlisted) > /dev/null 2>&1
			(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta_NODES_first_attmpt_shortlisted) > /dev/null 2>&1
			(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta) > /dev/null 2>&1
			(mkdir results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/refseq) > /dev/null 2>&1
			##---------------------------
			## coverage calculations
				grep ">" results/04_assembly/raw_files/$F1/$F1.fasta > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.strain_genome_coverage-calculation.tmp
				awk -F '_' '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.strain_genome_coverage-calculation.tmp > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.strain_genome_coverage.tmp
				rm results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.strain_genome_coverage-calculation.tmp 
			##---------------------------
			## get WGS-Plasmid-replicon
			echo "getting WGS-Plasmid-replicon for $F1 draft genome"
			blastn -db /media/swapnil/share/databases/plasmid_06022020_replicon/plasmid.fsa -query results/04_assembly/raw_files/$F1/$F1.fasta -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk '$14 > 90' > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.plasmid.blast.tmp
			WGS_pReplicon=$(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.plasmid.blast.tmp | awk -F'_' '{print $1}' | sort -u | tr "\n" "_" | sed '$s/.$//')
			if [ -z "$WGS_pReplicon" ]; then 
				WGS_pReplicon=$(echo "NA")
				else
				cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.plasmid.blast.tmp | awk '{print $2}' | sed -e "s/$/\.fa/" > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.WGS_pReplicon.contigs.list
			fi
			echo $WGS_pReplicon > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.WGS_pReplicon.list
			##---------------------------
			## get WGS-MLST-replicon
			echo "getting WGS-MLST-replicon for $F1 draft genome"
			blastn -db /media/swapnil/share/databases/plasmid_31052010_MLST_replicons/310519_all.fsa -query results/04_assembly/raw_files/$F1/$F1.fasta -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk '$14 > 90' > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.blast.tmp 
			WGS_replicon=$(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.blast.tmp | awk -F'_' '{print $1}' | sort -u | tr "\n" "_" | sed '$s/.$//')
			if [ -z "$WGS_replicon" ]; then 
				WGS_replicon=$(echo "NA")
				else
				cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.blast.tmp | awk '{print $2}' | sed -e "s/$/\.fa/" > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.WGS_replicon.contigs.list
			fi
			echo $WGS_replicon > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.WGS_replicon.list
			##---------------------------
			echo $F1 > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/strain_name.txt

			## blast contigs against database
			echo "BLAST $F1 contigs against plasmid_17112019_NCBI database"
			blastn -db /media/swapnil/share/databases/plasmid_17112019_NCBI/all_plasmid.fna -query results/04_assembly/raw_files/$F1/$F1.fasta -out results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.tmp -max_target_seqs 5 -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"
			awk '$15 > 25' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.tmp | awk '$12 < 300' > tmp && mv tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.tmp ;

			if grep -Fxq "sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.tmp; then
				:
				else
				ex -sc '1i|sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs' -cx results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.tmp
			fi
			awk '{print $1}' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.tmp | sed "s/ref|//g" |  sed "s/|//g" | sed '/sseqid/d' | sort -u > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.1_attempt_shortlisted.tmp
			awk '{print $2}' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.tmp | sed '/qseqid/d' | sort -u > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.nodes.1_attempt_shortlisted.tmp
			awk -f /home/swapnil/pipeline/tools/vlookup-plasmid-refseq.awk /media/swapnil/share/databases/plasmid_17112019_NCBI/all_plasmid.fna.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.1_attempt_shortlisted.tmp > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.annotation.tmp
			paste results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.1_attempt_shortlisted.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.annotation.tmp > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmids.csv

			#-------------------------------------------
			for ref_p in $(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.1_attempt_shortlisted.tmp);do
				echo "BLAST $F1 contigs against short-listed $ref_p plasmid"
				echo $ref_p > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_accession_numbers_unique.tmp.tmp
				awk '$1 == "'$ref_p'" {print $2}' /media/swapnil/share/databases/plasmid_17112019_NCBI/all_plasmid.fna.list > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_accession_numbers_unique.name.tmp

				perl /home/swapnil/pipeline/tools/fastagrep.pl -f results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_accession_numbers_unique.tmp.tmp /media/swapnil/share/databases/plasmid_17112019_NCBI/all_plasmid.fna > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/refseq/$ref_p.RefSeq.fasta

				grep $ref_p results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.tmp | awk '{print $2}' > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.tmp
				wc -l results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.tmp | awk '{print $1}' > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.count.tmp

				awk -F '_' '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.tmp > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.coverage.tmp
				
				perl /home/swapnil/pipeline/tools/fastagrep.pl -f results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.tmp results/04_assembly/raw_files/$F1/$F1.fasta > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/$F1.$ref_p.fasta

				cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.tmp | tr "\n" " " > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.2.tmp
				#---------------------

				sed 's/>.*/NNNNNNNNNN/g' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/$F1.$ref_p.fasta | sed "1i "'>'$F1.$ref_p"" > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/$F1.$ref_p.fasta.tmp 
				length=$(bioawk -c fastx '{ print $name, length($seq) }' <results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/$F1.$ref_p.fasta.tmp  | awk 'NR==1 {print $2}')

				if [ "$length" -gt "5000" ] ; then

					##---------------------------
					echo "getting pReplicon for $F1.$ref_p "
					blastn -db /media/swapnil/share/databases/plasmid_06022020_replicon/plasmid.fsa -query results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/$F1.$ref_p.fasta.tmp -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk '$14 > 90' > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.F_plasmid.blast.tmp 
					WGS_pReplicon=$(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.F_plasmid.blast.tmp | awk -F'_' '{print $1}' | sort -u | tr "\n" "_" | sed '$s/.$//')
					if [ -z "$WGS_pReplicon" ]; then 
						WGS_pReplicon=$(echo "NA")
						else
						cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.F_plasmid.blast.tmp | awk '{print $2}' | sed -e "s/$/\.fa/" > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.F_WGS_pReplicon.contigs.list
					fi
					echo $WGS_pReplicon > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.F_pReplicon.list

					##---------------------------
					echo "getting pMLST-replicon for $F1.$ref_p"
					blastn -db /media/swapnil/share/databases/plasmid_31052010_MLST_replicons/310519_all.fsa -query results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/$F1.$ref_p.fasta.tmp -max_hsps 1 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" | awk '$14 > 90' > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.F_blast.tmp
					WGS_replicon=$(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.F_blast.tmp | awk -F'_' '{print $1}' | sort -u | tr "\n" "_" | sed '$s/.$//')
					if [ -z "$WGS_replicon" ]; then 
						WGS_replicon=$(echo "NA")
						else
						cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.F_blast.tmp | awk '{print $2}' | sed -e "s/$/\.fa/" > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.F_WGS_replicon.contigs.list
					fi
					echo $WGS_replicon > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.F_replicon.list

					# calculate ANI ---------------------------------------------------------------

					mkdir /home/swapnil/tools/output
					( perl /home/swapnil/pipeline/tools/ANI_and_Total-aligned_plasmid_200.pl \
						--fd /home/swapnil/tools/blast-2.2.9-amd64-linux/bin/formatdb \
						-bl /home/swapnil/tools/blast-2.2.9-amd64-linux/bin/blastall \
						--qr results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/$F1.$ref_p.fasta.tmp \
						--sb results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/refseq/$ref_p.RefSeq.fasta \
						--od /home/swapnil/pipeline/tools/output >> results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp ) > /dev/null 2>&1

					sed -i "s/\.fasta\.tmp//g" results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp
					rm -r /home/swapnil/tools/output

					#------------------------------------------------------------------------------
					ANI_Cov=$(awk '{print $3}' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp)

					awk -F'/' '{print $NF}' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp > tmp && mv tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp ;
					sed -i 's/\.aligned\.joined\.fasta//g' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp ;

					awk '{queryCov= 100*$3 / $4} END {print queryCov}' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.query-covered.tmp ;
					awk '{refCov= 100*$3 / $5} END {print refCov}' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.ref-covered.tmp ;

					query_covered=$(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.query-covered.tmp | awk -F'.' '{print $1}' )
					refCov=$(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.ref-covered.tmp | awk -F'.' '{print $1}' )

					## if stat are enough good then carry out further analysis 
					if [ "$ANI_Cov" -gt "5000" ] && [ "$query_covered" -gt "25" ] && [ "$refCov" -gt "10" ]; then

						paste results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/strain_name.txt results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_accession_numbers_unique.tmp.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.query-covered.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.ref-covered.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.count.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.strain_genome_coverage.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.coverage.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.F_pReplicon.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.F_replicon.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.WGS_pReplicon.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.WGS_replicon.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_accession_numbers_unique.name.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.2.tmp >> results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_results.csv

						paste results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/strain_name.txt results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_accession_numbers_unique.tmp.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.query-covered.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.ANI.ref-covered.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.count.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.strain_genome_coverage.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.coverage.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.F_pReplicon.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.F_replicon.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.WGS_pReplicon.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.WGS_replicon.list results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_accession_numbers_unique.name.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.plasmid_accession_numbers_unique.nodes.2.tmp >> results/25_plasmid_finder_complete_BLAST"$s"/all.csv

						ssconvert results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.plasmid_results.csv results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/$F1.plasmid_results.csv.xlsx
						
						## --------------------------------------------------
						##  align contigs against Ref

							echo "aligning $F1 contigs against $ref_p"

							cp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/refseq/$ref_p.RefSeq.fasta /home/swapnil/tools/mauve/
							cp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/p_fasta/$F1.$ref_p.fasta /home/swapnil/tools/mauve/

							WorDir=$(echo $PWD)
							cd /home/swapnil/tools/mauve
							(java -Xmx50000m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output "$F1.$ref_p"_mauve -ref $ref_p.RefSeq.fasta  -draft $F1.$ref_p.fasta) > /dev/null 2>&1
							rm $ref_p.RefSeq.fasta
							rm $F1.$ref_p.fasta
							cd $WorDir

							mv /home/swapnil/tools/mauve/"$F1.$ref_p"_mauve results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/
							Var1=$(ls results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/"$F1.$ref_p"_mauve | sort -r | head -1)
							cp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/"$F1.$ref_p"_mauve/$Var1/$F1.$ref_p.fasta  results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.aligned.fasta
							cp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.aligned.fasta results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.aligned.joined.fasta
							sed -i 's/>.*/NNNNNNNNNN/g' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.aligned.joined.fasta
							sed -i "1i "'>'$F1.$ref_p"" results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.aligned.joined.fasta
							rm -r results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/"$F1.$ref_p"_mauve


						# Annotation of plasmid -------------------------------------------------------

						#if [ ! -d results/08_annotation_plasmid ]; then
						#mkdir results/08_annotation_plasmid
						#fi

						#if [ ! -d results/08_annotation_plasmid/$F1 ]; then
						#mkdir results/08_annotation_plasmid/$F1
						#fi

						#/home/swapnil/tools/prokka-1.13/bin/prokka --quiet --outdir results/08_annotation_plasmid/$F1/$F1.$ref_p --force --prefix $F1.$ref_p --addgenes --locustag $F1.$ref_p --strain $F1.$ref_p --rnammer results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.$ref_p.aligned.joined.fasta

						#------------------------------------------------------------------------------

					fi
				fi
			done
		done

		sed -i '1 i\strain plasmid strain-plasmid ANI Aprox_alignment_length query_total_length subject_total_length %-Query-Covered %-Ref-covered number_of_contigs Avg-Genome-Cov Avg-plasmid-Cov pReplicon pMLST WGS-pReplicon WGS-pMLST refSeq_Plasmid-name contigs' results/25_plasmid_finder_complete_BLAST"$s"/all.csv
		sed -i "s/ /\t/g" results/25_plasmid_finder_complete_BLAST"$s"/all.csv
		ssconvert results/25_plasmid_finder_complete_BLAST"$s"/all.csv results/25_plasmid_finder_complete_BLAST"$s"/all.csv.xlsx

		###############################################################################

		#echo "code to define and copy HC-plasmids are missing for below script"

		#cd results/25_plasmid_finder_complete_BLAST"$s"/all_HC_plasmid_fasta/
		#ls *.aligned.joined.fasta | sed 's/\.aligned\.joined\.fasta//g' > all_HC_plasmid_fasta.list.txt 
		#cd ..
		#cd ..
		#cd ..

		#mv results/25_plasmid_finder_complete_BLAST"$s"/all_HC_plasmid_fasta/all_HC_plasmid_fasta.list.txt results/25_plasmid_finder_complete_BLAST"$s"/

		#if [ ! -d results/25_plasmid_finder_complete_BLAST"$s"/replicon ]; then
		#mkdir results/25_plasmid_finder_complete_BLAST"$s"/replicon
		#fi

		#for plasmid in $(cat results/25_plasmid_finder_complete_BLAST"$s"/all_HC_plasmid_fasta.list.txt); do

			#echo "runnning... replicon search for $plasmid"

			#V1=$(echo $plasmid | awk -F'.' '{print $1}')


			#blastn -db /media/swapnil/share/databases/plasmid/plasmids.fasta -query results/08_annotation_plasmid/$V1/$plasmid/$plasmid.ffn -out results/25_plasmid_finder_complete_BLAST"$s"/replicon/$plasmid.replicon.tmp -evalue 1e-100 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs"

			#V2=$(awk -F'_' '{print $1}' results/25_plasmid_finder_complete_BLAST"$s"/replicon/$plasmid.replicon.tmp | uniq );

			#echo $plasmid  $V2 >> results/25_plasmid_finder_complete_BLAST"$s"/replicon/plasmid.replicons.csv

			#echo "finished... replicon search for $plasmid"
		#done


		## get shorlisted plasmids to a folder

		for shortlisted_p in $(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.1_attempt_shortlisted.tmp); do
			cp /media/swapnil/share/databases/plasmid_3/$shortlisted_p.fsa results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta_RefSeq_plasmid_first_attmpt_shortlisted/
		done

		#-------------------------------------------
		## get shortlisted nodes to a folder

		cp results/04_assembly/raw_files/$F1/$F1.fasta results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta/
		perl /home/swapnil/pipeline/tools/split_multifasta.pl -i results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta/$F1.fasta -o results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta/

		## from WGS-Plasmid-replicon
		awk '{print $2}' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.plasmid.blast.tmp > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.plasmid.blast.2.tmp
		## from WGS-MLST-replicon
		awk '{print $2}' results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.blast.tmp > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.blast.2.tmp

		cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.plasmid.blast.2.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/WGS/$F1.blast.2.tmp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.nodes.1_attempt_shortlisted.tmp | sort -u > results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.1_attempt_shortlisted.3.tmp

		for shortlisted_nodes in $(cat results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/$F1.RefseqPlasmid.1_attempt_shortlisted.3.tmp); do
			cp results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta/$shortlisted_nodes.fsa results/25_plasmid_finder_complete_BLAST"$s"/raw_files/$F1/tmp/fasta_NODES_first_attmpt_shortlisted/
		done
		echo "finished... Plasmid search for ILCC004 .........................."
	done


###############################################################################
echo "Completed.. plasmid serach complete ------------------------------------"
###############################################################################
