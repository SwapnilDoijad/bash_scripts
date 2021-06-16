#!/bin/bash
###############################################################################
#91 edirect 
###############################################################################
echo "started...... step-91 NCBI edirect utility  ----------------------------"
###############################################################################

PATH=$PATH:$HOME/tools/edirect

for F1 in $(cat list.assembly-acessions.txt); do
esearch -db nucleotide -query "$F1" | efetch -format docsum | grep "BioSample" | sed 's/<BioSample>//g' | sed 's/<\/BioSample>//g' | sed 's/	//g' >> list.biosample.txt
done 

(rm data/ncbi/all.tab) > /dev/null 2>&1
for F1 in $(cat list.biosample.txt); do

esearch -db biosample -query "$F1" | efetch -format native > $F1.tmp.tab

V1=$(grep "strain=" $F1.tmp.tab | sed "s/\/strain=//g" | sed "s/\"//g")
V2=$(grep "collection date=" $F1.tmp.tab | sed "s/\/collection date=//g" | sed "s/\"//g")
V3=$(grep "geographic location=" $F1.tmp.tab | sed "s/\/geographic location=//g" | sed "s/\"//g")
V4=$(grep "isolation source=" $F1.tmp.tab | sed "s/\/isolation source=//g" | sed "s/\"//g")
V5=$(grep "host=" $F1.tmp.tab | sed "s/\/host=//g" | sed "s/\"//g")
V6=$(grep "host disease=" $F1.tmp.tab | sed "s/\/host disease=//g" | sed "s/\"//g")
V7=$(grep "host age=" $F1.tmp.tab | sed "s/\/host age=//g" | sed "s/\"//g")
V8=$(grep "collected by=" $F1.tmp.tab | sed "s/\/collected by=//g" | sed "s/\"//g")
V9=$(grep "serovar=" $F1.tmp.tab | sed "s/\/serovar=//g" | sed "s/\"//g")
echo -e "$V1\t$V2\t$V3\t$V4\t$V5\t$V6\t$V7\t$V8\t$V9" >> all.tmp
sed -i "s/    //g" all.tmp
paste list.txt list.biosample.txt all.tmp
done

sed  -i '1i accession	biosample	strain	collection_date	geographic_location	isolation_source	host	host_disease	host_age	collected_by	serovar' all.tab

###############################################################################

echo "completed ...... step-91 NCBI e-direct utility -------------------------"

###############################################################################


