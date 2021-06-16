#!/bin/bash
###############################################################################
## download_assemblies_and_arrange_in_04_assembly
###############################################################################
## initial input and directory creation
    echo “Update: What is the name of the genus species?”
    read name
	echo "want to store the data on external drive, then provide full path? (for e.g. /media/swapnil/network/reads_database/p_Enterobacter_NCBI)"
	read path_tmp
	if [ ! -z $path_tmp ] ; then
		path=$(echo "$path_tmp/")
		else
        path=$(echo "$path_tmp")
	fi
	T=`date '+%d%m%Y_%H%M%S'`
	(mkdir data/ncbi)> /dev/null 2>&1
	(mkdir data/ncbi/tmp)> /dev/null 2>&1
	(mkdir data/ncbi/refseq)> /dev/null 2>&1
	(mkdir results)> /dev/null 2>&1
	(mkdir results/04_assembly)> /dev/null 2>&1
	(mkdir results/04_assembly/raw_files/)> /dev/null 2>&1
	(mkdir results/04_assembly/all_fasta)> /dev/null 2>&1
###############################################################################
## download NEW assemblies
## Download assemblies
for AssemblyBioSampleAccn in $(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt); do
	if [ -d "$path"data/ncbi/refseq/$AssemblyBioSampleAccn ] ; then
	echo "$AssemblyBioSampleAccn already downloaded"
	else
		echo "Downloading $AssemblyBioSampleAccn"
		## if RefSeq not available then loo
			downl_link_tmp1=$(awk -F'\t' '{print $23}' data/ncbi/metadata/tmp/$AssemblyBioSampleAccn.Assembly-BioSampleAccn.2.tmp)
		if [ "$downl_link_tmp1" == "NA" ]; then
			downl_link_tmp1=$(awk -F'\t' '{print $24}' data/ncbi/metadata/tmp/$AssemblyBioSampleAccn.Assembly-BioSampleAccn.2.tmp)
		fi 
		downl_link_ext=$(echo $downl_link_tmp1 | awk -F'/' '{print $NF}')
		downl_link=$(echo "$downl_link_tmp1"/"$downl_link_ext""_genomic.fna.gz")
		( mkdir data/ncbi/refseq/$AssemblyBioSampleAccn)> /dev/null 2>&1
		( wget $downl_link -P data/ncbi/refseq/$AssemblyBioSampleAccn/ ) > /dev/null 2>&1
		( rm data/ncbi/refseq/$AssemblyBioSampleAccn/*.fna ) > /dev/null 2>&1
		gunzip -f -k data/ncbi/refseq/$AssemblyBioSampleAccn/*.gz

		cp data/ncbi/refseq/$AssemblyBioSampleAccn/*.fna data/ncbi/refseq/$AssemblyBioSampleAccn/$AssemblyBioSampleAccn.fna
		cp data/ncbi/refseq/$AssemblyBioSampleAccn/$AssemblyBioSampleAccn.fna data/ncbi/refseq/$AssemblyBioSampleAccn/$AssemblyBioSampleAccn.joined.fna
		sed -i 's/>.*/NNNNNNNNNN/g' data/ncbi/refseq/$AssemblyBioSampleAccn/$AssemblyBioSampleAccn.joined.fna
		sed -i "1i "'>'$AssemblyBioSampleAccn"" data/ncbi/refseq/$AssemblyBioSampleAccn/$AssemblyBioSampleAccn.joined.fna
		(mkdir results/04_assembly/raw_files/$AssemblyBioSampleAccn ) > /dev/null 2>&1
		cp data/ncbi/refseq/$AssemblyBioSampleAccn/$AssemblyBioSampleAccn.joined.fna results/04_assembly/raw_files/$AssemblyBioSampleAccn/$AssemblyBioSampleAccn.joined.fasta
		cp results/04_assembly/raw_files/$AssemblyBioSampleAccn/$AssemblyBioSampleAccn.joined.fasta results/04_assembly/all_fasta/$AssemblyBioSampleAccn.fasta
		if [ "$(ls -A results/04_assembly/raw_files/$AssemblyBioSampleAccn)" ]; then
		:
		else
		rm -rf results/04_assembly/raw_files/$AssemblyBioSampleAccn 
		fi

		if [ "$(ls -A "$path"data/ncbi/refseq/$AssemblyBioSampleAccn )" ]; then
		:
		else
		rm -rf results/04_assembly/raw_files/$AssemblyBioSampleAccn 
		fi
		## to move the final files to network/external drive 
		if [ ! -z $path_tmp ] ; then
			echo moving files to $path
			(mkdir $path/data)> /dev/null 2>&1
			(mkdir $path/data/ncbi)> /dev/null 2>&1
			(mkdir $path/data/ncbi/refseq)> /dev/null 2>&1
			(mkdir $path/data/ncbi/refseq/$AssemblyBioSampleAccn)> /dev/null 2>&1
			cp -rf data/ncbi/refseq/$AssemblyBioSampleAccn $path/data/ncbi/refseq/
			rm -rf data/ncbi/refseq/$AssemblyBioSampleAccn
		fi
	fi
done
exit
###############################################################################
## map contigs against reference

#echo "Do you want to map the contigs against the reference? yes or PRESS enter to skip the mapping"
#read answer

#if [ "$answer" == "yes" ] ; then
#cp results/00_ref/ref.*.fasta /home/swapnil/tools/mauve/

#for F1 in $(cat data/ncbi/list.assemblies.accessions.downloaded.tmp)
#do
#echo "running.... step-6 mapping filtered contigs to the reference for..." $F1

#cp results/04_assembly/raw_files/$F1/contigs.fasta /home/swapnil/tools/mauve/$F1.contigs.fasta

#WorDir=$(echo $PWD)

#cd /home/swapnil/tools/mauve
#(java -Xmx5000m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output $F1 -ref ref.*.fasta -draft $F1.contigs.fasta) > /dev/null 2>&1 

#cd $WorDir
#Var1=$(ls /home/swapnil/tools/mauve/$F1/ | sort -r | head -1)
#cp /home/swapnil/tools/mauve/$F1/$Var1/$F1.contigs.fasta  results/04_assembly/raw_files/$F1/"$F1".aligned.fasta
#cp results/04_assembly/raw_files/$F1/"$F1".aligned.fasta results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta
#sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta
#sed -i "1i "'>'$F1"" results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta

#cp results/04_assembly/raw_files/$F1/"$F1".aligned.joined.fasta results/04_assembly/all_fasta/$F1.fasta

#rm -rf /home/swapnil/tools/mauve/$F1
#rm /home/swapnil/tools/mauve/$F1.contigs.fasta

#done
#fi
#------------------------------------------------------------------------------
