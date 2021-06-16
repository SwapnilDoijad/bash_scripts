#!/bin/bash
###############################################################################
#91 edirect 
## https://www.ncbi.nlm.nih.gov/books/NBK179288/
###############################################################################

echo "started...... step-91 NCBI edirect utility  ----------------------------"

###############################################################################

PATH=$PATH:$HOME/tools/edirect

if [ ! -d results/91_NCBI_metadata_of_biosample ]; then
mkdir results/91_NCBI_metadata_of_biosample
fi

for F1 in $(cat list.SRA.txt); do

## understand the four level hierarchial structure of output docusm or xml and choose accordingly
# Exploration Argument Hierarchy
#    -pattern         Name of record within set
#    -group             Use of different argument
#    -block               names allows command-line
#    -subset                control of nested looping
# esearch -db biosample -query "Listeria ivanovii" | efetch -format docsum  | xtract -pattern DocumentSummary -group SampleData -block Comment -element Paragra
## this will give multiple attributes within <> </>
# esearch -db biosample -query "Listeria ivanovii" | efetch -format docsum  | xtract -pattern DocumentSummary -group SampleData -block Attributes -element Attribute
## if comment
# esearch -db biosample -query "Listeria ivanovii" | efetch -format docsum  | xtract -pattern DocumentSummary -group SampleData -block Organism -if @taxonomy_id -equals "202751" -element OrganismName
# esearch -db biosample -query "Listeria ivanovii" | efetch -format docsum  | xtract -pattern DocumentSummary -group SampleData -block BioSample -if @publication_date -equals "2014-08-06T00:00:00.000" -subset Ids -element Id
# esearch -db biosample -query "Listeria ivanovii" | efetch -format docsum  | xtract -pattern DocumentSummary -element SourceSample -group SampleData -block BioSample -if @publication_date -equals "2013-02-12T00:00:00.000" -subset Attributes -element Attribute
## countrywise ##filter
# echo $( esearch -db biosample -query "Listeria monocytogenes" | efilter -query "India" | efetch -format docsum | xtract -pattern DocumentSummary -element SourceSample _) | sed 's/BioSample:/\n/g' | sed '/^$/d'
## ##filter 'and' condition
# echo $( esearch -db biosample -query "Listeria monocytogenes" | efilter -query "India AND Buttermilk" | efetch -format docsum | xtract -pattern DocumentSummary -element SourceSample _) | sed 's/BioSample:/\n/g' | sed '/^$/d'

# elink
# esearch -db biosample -query "SAMN06109051" | elink -target assembly | efetch -format docsum 

# Download cds or fasta
# esearch -db nucleotide -query NC_003210 | efetch -format fasta_cds_aa

# get biosample from NCBI genome accession
# esearch -db nucleotide -query NZ_CP018926 | efetch -format docsum | grep "BioSample" | sed 's/<BioSample>//g' | sed 's/<\/BioSample>//g' | sed 's/	//g'
# esearch -db biosample -query "Listeria ivanovii" | efetch -format docsum  | xtract -pattern DocumentSummary -element SourceSample

# get pubmed data
# esearch -db pubmed -query "lycopene cyclase" | efetch -format abstract > results/91_NCBI_metadata_of_biosample/$F1.pubmed.tab

# bioproject data
# esearch -db bioproject -query "$F1" | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -block Attribute -element Attribute > results/91_NCBI_metadata_of_biosample/$F1.bioproject.biosample.tab

# esearch -db biosample -query "$F1" | efetch -format docsum | xtract -pattern DocumentSummary -block Attribute -element Attribute >> results/91_NCBI_metadata_of_biosample/biosample.tab

# esearch -db biosample -query "$F1" | efetch -format native > results/91_NCBI_metadata_of_biosample/$F1.tmp.tab

# How many Lm from India ? 
# esearch -db biosample -query "Listeria monocytogenes" | efetch -format native | grep -c "India"

# Download sequences from NCBI using edirect using bioproject accession or ID
# esearch -db bioproject -query "PRJNA285593" | elink -target nuccore | efetch -format fasta

# Get all CDS from a genome
# esearch -db protein -query 302315370| elink -target nuccore |efetch -format ft| grep -A 4 --no-group-separator CDS

## extracting Abr table  from xml to csv
# esearch -db biosample -query "SAMN05201844" | efetch -format docsum | xtract -pattern DocumentSummary -group Comment -block Row -sep "|" -element Cell | sed 's/\t/\n/g' | sed 's/|/\t/g'

## input SRA, links to biosample and get all possible meadata (then needs to define what needs to be extracted for e.g strain, host, date etc.)
esearch -db SRA -query "$F1" | elink -target biosample | efetch -format native > results/91_NCBI_metadata_of_biosample/$F1.tmp.tab

V1=$(grep "strain=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/strain=//g" | sed "s/\"//g")
V2=$(grep "collection date=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/collection date=//g" | sed "s/\"//g")
V3=$(grep "geographic location=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/geographic location=//g" | sed "s/\"//g")
V4=$(grep "isolation source=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/isolation source=//g" | sed "s/\"//g")
V5=$(grep "host=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/host=//g" | sed "s/\"//g")
V6=$(grep "host disease=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/host disease=//g" | sed "s/\"//g")
V7=$(grep "host age=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/host age=//g" | sed "s/\"//g")
V8=$(grep "collected by=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/collected by=//g" | sed "s/\"//g")
V9=$(grep "serovar=" results/91_NCBI_metadata_of_biosample/$F1.tmp.tab | sed "s/\/serovar=//g" | sed "s/\"//g")
echo -e "$V1\t$V2\t$V3\t$V4\t$V5\t$V6\t$V7\t$V8\t$V9" >> results/91_NCBI_metadata_of_biosample/all.tab
sed -i "s/    //g" results/91_NCBI_metadata_of_biosample/all.tab
done

sed  -i '1i strain	collection_date	geographic_location	isolation_source	host	host_disease	host_age	collected_by	serovar' results/91_NCBI_metadata_of_biosample/all.tab

###############################################################################

echo "completed ...... step-91 NCBI e-direct utility -------------------------"

###############################################################################


