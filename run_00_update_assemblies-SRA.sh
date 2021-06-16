#!/bin/bash
###############################################################################
PATH=$PATH:$HOME/tools/edirect

if [ ! -d data/update ]; then
mkdir data/update
fi

if [ ! -d data/update/tmp ]; then
mkdir data/update/tmp
fi

#------------------------------------------------------------------------------
echo “Update: What is the name of the genus species?”
read name

echo "Checking new assemblies for $name"
(esearch -db assembly -query "$name" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession) > data/tmp/"$name".AssemblyAccession.to-update.tmp
cat data/tmp/"$name".AssemblyAccession.tmp data/tmp/"$name".AssemblyAccession.to-update.tmp | sed '/^\s*$/d' | sort | uniq -u > data/list."$name".AssemblyAccession.to-update.txt
V1=$(wc -l data/list."$name".AssemblyAccession.to-update.txt | awk '{print $1}')
echo "$V1 new assemblies are available, respective new accession numbers are placed in in data/list."$name".AssemblyAccession.to-update.txt"

echo "Checking new SRA for $name"
(esearch -db sra -query "$name" | esummary -format runinfo -mode xml | xtract -pattern Row -element Run) > data/tmp/"$name".SRA-run.to-update.tmp
cat data/tmp/"$name".SRA-run.tmp data/tmp/"$name".SRA-run.to-update.tmp | sed '/^\s*$/d' | sort | uniq -u > data/list."$name".SRA-run.to-update.txt
V2=$(wc -l data/list."$name".SRA-run.to-update.txt | awk '{print $1}')
echo "$V2 new SRA are available, respective new accession numbers are placed in in data/list."$name".SRA-run.to-update.txt"

###############################################################################

echo "Do you want to get the metadata for new assemblies and SRA? answer yes or PRESS ENTER to skip" 
read F2
if [ "$F2" == "yes" ]; then

#------------------------------------------------------------------------------
for F1 in $(cat data/list."$name".AssemblyAccession.to-update.txt); do
echo "extracting metadata for $F1"
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,AssemblyAccession) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".AssemblyAccession.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Organism) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".Organism.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,SpeciesName) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".SpeciesName.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Sub_type) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".Sub_type.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Sub_value) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".Sub_value.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Isolate) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".Isolate.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,AssemblyName) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".AssemblyName.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,AssemblyStatus) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".AssemblyStatus.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,WGS) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".WGS.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Coverage) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".Coverage.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,BioSampleAccn) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".BioSampleAccn.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,SubmitterOrganization) | awk -F'\t' '{print $2}' >data/update/tmp/"$F1".SubmitterOrganization.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Taxid) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".Taxid.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,SubmissionDate) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".SubmissionDate.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,AsmReleaseDate_GenBank) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".AsmReleaseDate_GenBank.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,Date_GenBank) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".Date_GenBank.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,LastUpdateDate) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".LastUpdateDate.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,FtpPath_RefSeq) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".FtpPath_RefSeq.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession,FtpPath_GenBank) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".FtpPath_GenBank.tmp

(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession -block Stat -if "@category" -equals total_length -element Stat) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".total_length.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession -block Stat -if "@category" -equals contig_count -element Stat) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".contig_count.tmp
(esearch -db assembly -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -element  AssemblyAccession -block Stat -if "@category" -equals contig_n50 -element Stat) | awk -F'\t' '{print $2}' > data/update/tmp/"$F1".contig_n50.tmp

cd data/update/tmp/
paste "$F1".AssemblyAccession.tmp "$F1".total_length.tmp "$F1".contig_count.tmp "$F1".contig_n50.tmp "$F1".Organism.tmp "$F1".SpeciesName.tmp "$F1".Sub_type.tmp "$F1".Sub_value.tmp "$F1".Isolate.tmp "$F1".AssemblyName.tmp "$F1".AssemblyStatus.tmp "$F1".WGS.tmp "$F1".Coverage.tmp "$F1".BioSampleAccn.tmp "$F1".SubmitterOrganization.tmp "$F1".Taxid.tmp "$F1".SubmissionDate.tmp "$F1".AsmReleaseDate_GenBank.tmp "$F1".LastUpdateDate.tmp "$F1".FtpPath_RefSeq.tmp "$F1".FtpPath_GenBank.tmp > 00_"$F1"_assemblies.metadata.slow.tab
cd ..
cd ..
cd ..

done

cat data/update/tmp/*assemblies.metadata.slow.tab > data/update/update.assemblies.metadata.slow.tab
sed  -i '1i AssemblyAccession	total_length	contig_count	contig_n50	Organism	SpeciesName	Sub_type	Sub_value	Isolate	AssemblyName	AssemblyStatus	WGS	Coverage	BioSampleAccn	SubmitterOrganization	Taxid	SubmissionDate	AsmReleaseDate_GenBank	LastUpdateDate	FtpPath_RefSeq	FtpPath_GenBank' data/update/update.assemblies.metadata.slow.tab
###############################################################################
## SRA-metadata 

for F1 in $(cat data/list."$name".SRA-run.to-update.txt); do
(esearch -db sra -query "$F1" | esummary -format runinfo -mode xml | xtract -pattern Row -element Run,bases,avgLength,size_MB,Experiment,LibrarySource,LibraryLayout,Platform,Model,SRAStudy,BioProject,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,Submission,ReleaseDate,download_path) > data/update/tmp/00_"$F1"_SRA.available-to-download.tab

#------------------------------------------------------------------------------
## SRA--biosample-metadata: input SRA, links to biosample and get all possible meadata (then needs to define what needs to be extracted for e.g strain, host, date etc.)

echo "getting SRA's biosample metadata for $F1"

awk -F'\t' '{print $14}' data/update/tmp/00_"$F1"_SRA.available-to-download.tab > data/tmp/$F1.list.biosample.accessions.txt

for F3 in $(cat data/tmp/$F1.list.biosample.accessions.txt); do

esearch -db biosample -query "$F3" | efetch -format native > data/update/tmp/$F1.$F3.tmp.tab

V1=$(grep "strain=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/strain=//g" | sed "s/\"//g")
V2=$(grep "collection date=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/collection date=//g" | sed "s/\"//g")
V3=$(grep "geographic location=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/geographic location=//g" | sed "s/\"//g")
V4=$(grep "isolation source=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/isolation source=//g" | sed "s/\"//g")
V5=$(grep "host=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/host=//g" | sed "s/\"//g")
V6=$(grep "host disease=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/host disease=//g" | sed "s/\"//g")
V7=$(grep "host age=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/host age=//g" | sed "s/\"//g")
V8=$(grep "collected by=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/collected by=//g" | sed "s/\"//g")
V9=$(grep "serovar=" data/update/tmp/$F1.$F3.tmp.tab | sed "s/\/serovar=//g" | sed "s/\"//g")
V10=$(sed -n -e 's/^.*BioSample: //p' data/update/tmp/$F1.$F3.tmp.tab | sed 's/;//g' | awk '{print $1}')
V11=$(sed -n -e 's/^.*SRA: //p' data/update/tmp/$F1.$F3.tmp.tab | sed 's/;//g' | awk '{print $1}')

echo -e "$V11\t$V10\t$V1\t$V2\t$V3\t$V4\t$V5\t$V6\t$V7\t$V8\t$V9" > data/update/tmp/00_"$F1"_SRA.available-to-download.metadata.tab
sed -i "s/    //g" data/update/tmp/00_"$F1"_SRA.available-to-download.metadata.tab
done
done
#------------------------------------------------------------------------------
cat data/update/tmp/*SRA.available-to-download.tab > data/update/tmp/update.SRA.available-to-download.tab1
sed -i '1i	Run	bases	avgLength	size_MB	Experiment	LibrarySource	LibraryLayout	Platform	Model	SRAStudy	BioProject	ProjectID	Sample	BioSample	SampleType	TaxID	ScientificName	SampleName	Submission	ReleaseDate	download_path' data/update/tmp/update.SRA.available-to-download.tab1

cat data/update/tmp/*SRA.available-to-download.metadata.tab > data/update/tmp/update.SRA.available-to-download.metadata.tab1 ;
sed -i '1i	SRA	BioSample	strain	collection_date	geographic_location	isolation_source	host	host_disease	host_age	collected_by	serovar' data/update/tmp/update.SRA.available-to-download.metadata.tab1
#------------------------------------------------------------------------------

paste data/update/tmp/update.SRA.available-to-download.metadata.tab1 data/update/tmp/update.SRA.available-to-download.tab1 > data/update/update.SRA.metadata.tab


else
:
fi
