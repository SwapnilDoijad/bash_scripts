#!/bin/bash
###############################################################################
## get the assemblies/SRA metadata and accession to downloads 

###############################################################################

echo "started...... step-00 get the assemblies/SRA metadata and accession ...."

###############################################################################
PATH=$PATH:$HOME/tools/edirect

if [ ! -d data ]; then
mkdir data
fi

if [ ! -d data/tmp ]; then
mkdir data/tmp
fi


#------------------------------------------------------------------------------
## assemblies

echo “What is the name of the genus species?”
read F1

echo "Do you want to run slow and accurate process? answer yes or PRESS ENTER to skip to fast process" 
read F2
if [ "$F2" == "yes" ]; then

SECONDS=0
echo "OK, getting assembly accession and metadata for $F1 (Process: slow and accurate-output)"

( mkdir data/tmp/tmp-assembly ) > /dev/null 2>&1

(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,AssemblyAccession) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".AssemblyAccession.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Organism) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".Organism.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,SpeciesName) |  awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".SpeciesName.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Sub_type) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".Sub_type.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Sub_value) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".Sub_value.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Isolate) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".Isolate.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,AssemblyName) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".AssemblyName.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,AssemblyStatus) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".AssemblyStatus.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,WGS) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".WGS.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Coverage) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".Coverage.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,BioSampleAccn) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".BioSampleAccn.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,SubmitterOrganization) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".SubmitterOrganization.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Taxid) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".Taxid.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,SubmissionDate) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".SubmissionDate.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,AsmReleaseDate_GenBank) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".AsmReleaseDate_GenBank.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Date_GenBank) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".Date_GenBank.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,LastUpdateDate) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".LastUpdateDate.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,FtpPath_RefSeq) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".FtpPath_RefSeq.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,FtpPath_GenBank) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".FtpPath_GenBank.tmp

(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession -block Stat -if "@category" -equals total_length -element Stat) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".total_length.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession -block Stat -if "@category" -equals contig_count -element Stat) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".contig_count.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession -block Stat -if "@category" -equals contig_n50 -element Stat) | awk -F'\t' '{print $2}' > data/tmp/tmp-assembly/"$F1".contig_n50.tmp

cd data/tmp/tmp-assembly/
paste "$F1".AssemblyAccession.tmp "$F1".total_length.tmp "$F1".contig_count.tmp "$F1".contig_n50.tmp "$F1".Organism.tmp "$F1".SpeciesName.tmp "$F1".Sub_type.tmp "$F1".Sub_value.tmp "$F1".Isolate.tmp "$F1".AssemblyName.tmp "$F1".AssemblyStatus.tmp "$F1".WGS.tmp "$F1".Coverage.tmp "$F1".BioSampleAccn.tmp "$F1".SubmitterOrganization.tmp "$F1".Taxid.tmp "$F1".SubmissionDate.tmp "$F1".AsmReleaseDate_GenBank.tmp "$F1".LastUpdateDate.tmp "$F1".FtpPath_RefSeq.tmp "$F1".FtpPath_GenBank.tmp > 00_"$F1".assemblies.metadata.slow.tmp
cd ..
cd ..
cd ..

sed  -i '1i\AssemblyAccession	total_length	contig_count	contig_n50	Organism	SpeciesName	Sub_type	Sub_value	Isolate	AssemblyName	AssemblyStatus	WGS	Coverage	BioSampleAccn	SubmitterOrganization	Taxid	SubmissionDate	AsmReleaseDate_GenBank	LastUpdateDate	FtpPath_RefSeq	FtpPath_GenBank' data/tmp/tmp-assembly/00_"$F1".assemblies.metadata.slow.tmp


awk -F'\t' 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' data/tmp/tmp-assembly/00_"$F1".assemblies.metadata.slow.tmp > data/tmp/00_"$F1".assemblies.metadata.slow.tab
#------------------------------------------------------------------------------
##  Assembly biosample-metadata: input biosample and get all possible meadata 

echo "getting biosample-metadata for assemblies"

(mkdir data/tmp/tmp-biosample) > /dev/null 2>&1

(rm data/tmp/tmp-biosample/"$F1".list.assemblies.biosample.accessions.metadata.tmp) > /dev/null 2>&1

awk -F'\t' 'NR>1 {print $14}' data/tmp/00_"$F1".assemblies.metadata.slow.tab > data/tmp/tmp-biosample/"$F1".list.assemblies.biosample.accessions.tmp

for F3 in $(cat data/tmp/tmp-biosample/"$F1".list.assemblies.biosample.accessions.tmp); do

esearch -db biosample -query "$F3" | efetch -format native > data/tmp/tmp-biosample/"$F1".$F3.tmp.tab

V1=$(grep "strain=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/strain=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V2=$(grep "collection date=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/collection date=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V3=$(grep "geographic location=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/geographic location=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V4=$(grep "isolation source=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/isolation source=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V5=$(grep "host=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/host=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V6=$(grep "host disease=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/host disease=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V7=$(grep "host age=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/host age=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V8=$(grep "collected by=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/collected by=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V9=$(grep "serovar=" data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed "s/\/serovar=//g" | sed "s/\"//g" | awk -F'\t' 'NR==1{print $1}')
V10=$(sed -n -e 's/^.*BioSample: //p' data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed 's/;//g' | awk '{print $1}' | awk -F'\t' 'NR==1{print $1}')
V11=$(sed -n -e 's/^.*SRA: //p' data/tmp/tmp-biosample/"$F1".$F3.tmp.tab | sed 's/;//g' | awk '{print $1}' | awk -F'\t' 'NR==1{print $1}')

echo -e "$V11\t$V10\t$V1\t$V2\t$V3\t$V4\t$V5\t$V6\t$V7\t$V8\t$V9" >> data/tmp/tmp-biosample/"$F1".list.assemblies.biosample.accessions.metadata.tmp
done
sed -i "s/    //g" data/tmp/tmp-biosample/"$F1".list.assemblies.biosample.accessions.metadata.tmp

sed  -i '1i SRA	BioSample	strain	collection_date	geographic_location	isolation_source	host	host_disease	host_age	collected_by	serovar' data/tmp/tmp-biosample/"$F1".list.assemblies.biosample.accessions.metadata.tmp

awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' data/tmp/tmp-biosample/"$F1".list.assemblies.biosample.accessions.metadata.tmp > data/tmp/00_"$F1".list.assemblies.biosample.accessions.metadata.tab

#----------------------

paste data/tmp/00_"$F1".assemblies.metadata.slow.tab data/tmp/00_"$F1".list.assemblies.biosample.accessions.metadata.tab > data/00_"$F1".assemblies.metadata.tab


duration=$SECONDS
echo "finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds" 
exit
#------------------------------------------------------------------------------
## for fast processing


else

SECONDS=0

echo "OK, getting assembly accession and metadata for $F1 (Process: fast but check the output for column mismatch)"

esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element Organism,SpeciesName,Sub_type,Sub_value,Isolate,AssemblyName,AssemblyStatus,WGS,Coverage,BioSampleAccn,SubmitterOrganization,Taxid,SubmissionDate,AsmReleaseDate_GenBank,LastUpdateDate,FtpPath_RefSeq,FtpPath_GenBank > data/tmp/00_"$F1"_assemblies.available-to-download.fast.1.tmp

awk -F'\t' 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' data/tmp/00_"$F1"_assemblies.available-to-download.fast.1.tmp > data/tmp/00_"$F1"_assemblies.available-to-download.fast.2.tmp

esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession -block Stat -if "@category" -equals total_length -element Stat -block Stat -if "@category" -equals contig_count -element Stat -block Stat -if "@category" -equals contig_n50 -element Stat > data/tmp/00_"$F1"_assemblies.available-to-download.fast.3.tmp

awk -F'\t' 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' data/tmp/00_"$F1"_assemblies.available-to-download.fast.3.tmp > data/tmp/00_"$F1"_assemblies.available-to-download.fast.4.tmp

paste data/tmp/00_"$F1"_assemblies.available-to-download.fast.4.tmp data/tmp/00_"$F1"_assemblies.available-to-download.fast.2.tmp  > data/00_"$F1"_assemblies.metadata.fast.tab

sed  -i '1i AssemblyAccession	total_length	contig_count	contig_n50	Organism	SpeciesName	Sub_type	Sub_value	Isolate	AssemblyName	AssemblyStatus	WGS	Coverage	BioSampleAccn	SubmitterOrganization	Taxid	SubmissionDate	AsmReleaseDate_GenBank	LastUpdateDate	FtpPath_RefSeq	FtpPath_GenBank' data/00_"$F1"_assemblies.metadata.fast.tab

duration=$SECONDS
echo "finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds" 

fi
###############################################################################
## SRA-metadata 

echo "Should SRA metadata also be extracted for $F1 answer yes or PRESS ENTER to skip" 
read F2
if [ "$F2" == "yes" ]; then

echo "getting SRA's metadata for $F1"

esearch -db sra -query "$F1" | esummary -format runinfo -mode xml | xtract -pattern Row -element Run,bases,avgLength,size_MB,Experiment,LibrarySource,LibraryLayout,Platform,Model,SRAStudy,BioProject,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,Submission,ReleaseDate,download_path > data/tmp/00_"$F1"_SRA.available-to-download.tab

sed -i '1i	Run	bases	avgLength	size_MB	Experiment	LibrarySource	LibraryLayout	Platform	Model	SRAStudy	BioProject	ProjectID	Sample	BioSample	SampleType	TaxID	ScientificName	SampleName	Submission	ReleaseDate	download_path' data/tmp/00_"$F1"_SRA.available-to-download.tab

awk -F'\t' 'NR>1 {print $1}' data/tmp/00_"$F1"_SRA.available-to-download.tab > data/tmp/"$F1".SRA-run.tmp
#------------------------------------------------------------------------------
## SRA--biosample-metadata: input SRA, links to biosample and get all possible meadata (then needs to define what needs to be extracted for e.g strain, host, date etc.)

echo "getting SRA's biosample metadata for $F1"
(rm data/tmp/00_"$F1"_SRA.available-to-download.metadata.tab) > /dev/null 2>&1

awk -F'\t' 'NR>1 {print $14}' data/tmp/00_"$F1"_SRA.available-to-download.tab > data/tmp/$F1.list.biosample.accessions.txt

for F3 in $(cat data/tmp/$F1.list.biosample.accessions.txt); do

esearch -db biosample -query "$F3" | efetch -format native > data/tmp/$F1.$F3.tmp.tab

V1=$(grep "strain=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/strain=//g" | sed "s/\"//g")
V2=$(grep "collection date=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/collection date=//g" | sed "s/\"//g")
V3=$(grep "geographic location=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/geographic location=//g" | sed "s/\"//g")
V4=$(grep "isolation source=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/isolation source=//g" | sed "s/\"//g")
V5=$(grep "host=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/host=//g" | sed "s/\"//g")
V6=$(grep "host disease=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/host disease=//g" | sed "s/\"//g")
V7=$(grep "host age=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/host age=//g" | sed "s/\"//g")
V8=$(grep "collected by=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/collected by=//g" | sed "s/\"//g")
V9=$(grep "serovar=" data/tmp/$F1.$F3.tmp.tab | sed "s/\/serovar=//g" | sed "s/\"//g")
V10=$(sed -n -e 's/^.*BioSample: //p' data/tmp/$F1.$F3.tmp.tab | sed 's/;//g' | awk '{print $1}')
V11=$(sed -n -e 's/^.*SRA: //p' data/tmp/$F1.$F3.tmp.tab | sed 's/;//g' | awk '{print $1}')

echo -e "$V11\t$V10\t$V1\t$V2\t$V3\t$V4\t$V5\t$V6\t$V7\t$V8\t$V9" >> data/tmp/00_"$F1"_SRA.available-to-download.metadata.tab
sed -i "s/    //g" data/tmp/00_"$F1"_SRA.available-to-download.metadata.tab
done

sed  -i '1i SRA	BioSample	strain	collection_date	geographic_location	isolation_source	host	host_disease	host_age	collected_by	serovar' data/tmp/00_"$F1"_SRA.available-to-download.metadata.tab

paste data/tmp/00_"$F1"_SRA.available-to-download.tab data/tmp/00_"$F1"_SRA.available-to-download.metadata.tab > data/tmp/00_"$F1"_SRA.metadata.tab

awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' data/tmp/00_"$F1"_SRA.metadata.tab > data/00_"$F1"_SRA.metadata.tab

else
:
fi

###############################################################################

echo "completed ...... step-00 get the assemblies/SRA metadata and accession.."

###############################################################################




