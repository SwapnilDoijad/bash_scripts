#!/bin/bash
###############################################################################
#91 edirect 
## ncbi-genome-download -n -M type -g Klebsiella bacteria
###############################################################################
echo "started...... step-91 NCBI edirect utility  ----------------------------"
###############################################################################
## inital file preparation and directory creation
    T=`date '+%d%m%Y_%H%M%S'`

    #echo "Do you want to download these files and arrange in 04_assembly folder? type y or press ENTER"
    #read download_ans
    download_ans=y
    #if [ ! -z $download_ans ] ; then
    #echo "Are these type strains? type y or press ENTER"
    #read type_ans
    #fi

    echo “Update: What is the name of the genus species?”
    read name

    PATH=$PATH:$HOME/tools/edirect

    (mkdir data) > /dev/null 2>&1
    (mkdir data/ncbi) > /dev/null 2>&1
    (mkdir data/ncbi/tmp) > /dev/null 2>&1
    (mkdir data/ncbi/metadata) > /dev/null 2>&1 
    (mkdir data/ncbi/metadata/tmp) > /dev/null 2>&1 
    (touch data/ncbi/tmp/list.biosample.txt) > /dev/null 2>&1
    (touch data/ncbi/tmp/all.tmp) > /dev/null 2>&1
###############################################################################
## accession to biosample 
    ## check if these not present previously

    (touch data/ncbi/tmp/"$name".BioSampleAccn."$T".tmp ) > /dev/null 2>&1
    for accession in $(cat list.accession.txt); do
        bioSample=$( esearch -db nucleotide -query "$accession" | efetch -format docsum | xtract -pattern  DocumentSummarySet -element BioSample | head -1 )
        echo $bioSample
        ## check if biosample empty
        if [ ! -z $bioSample ] ; then
            if grep -q $bioSample data/ncbi/tmp/"$name".BioSampleAccn."$T".tmp ; then
                :
                else
                echo $bioSample >> data/ncbi/tmp/"$name".BioSampleAccn."$T".tmp
                echo $accession $bioSample >> data/ncbi/tmp/"$name".accession.BioSampleAccn."$T".tmp
            fi
            else
            echo $accession "NA" >> data/ncbi/tmp/"$name".accession.BioSampleAccn."$T".tmp
        fi
    done
    
    (touch data/ncbi/metadata/tmp/old_BioSampleAccn.txt)> /dev/null 2>&1
    (ls data/ncbi/metadata/tmp/*.Assembly-BioSampleAccn.tmp | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' | sort > data/ncbi/metadata/tmp/old_BioSampleAccn.txt) > /dev/null 2>&1
    (cat data/ncbi/tmp/"$name".BioSampleAccn."$T".tmp | sort -u > data/ncbi/metadata/tmp/new_BioSampleAccn.txt)> /dev/null 2>&1
    comm -13 data/ncbi/metadata/tmp/old_BioSampleAccn.txt data/ncbi/metadata/tmp/new_BioSampleAccn.txt > data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt
    cat data/ncbi/metadata/tmp/old_BioSampleAccn.txt data/ncbi/metadata/tmp/new_BioSampleAccn.txt | sort | uniq > data/ncbi/list."$name".Assembly-BioSampleAccn.all-BioSampleAccn.txt
###############################################################################
## metadata
    for AssemblyBiosampleAcc in $(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt); do
        echo "extracting Assembly-BioSampleAccn metadata for $AssemblyBiosampleAcc"
        esearch -db assembly -query "$AssemblyBiosampleAcc" | efetch -format docsum > data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp
        AssemblyAccession=$(grep AssemblyAccession data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 |  awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        Organism=$(grep Organism data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        SpeciesName=$(grep SpeciesName data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        Sub_type=$(grep Sub_type data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        Sub_value=$(grep Sub_value data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        Isolate=$(grep Isolate data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        AssemblyName=$(grep AssemblyName data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        AssemblyStatus=$(grep AssemblyStatus data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        WGS=$(grep WGS data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        Coverage=$(grep Coverage data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        SubmitterOrganization=$(grep SubmitterOrganization data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        Taxid=$(grep Taxid data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        SubmissionDate=$(grep SubmissionDate data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        AsmReleaseDate_GenBank=$(grep AsmReleaseDate_GenBank data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        Date_GenBank=$(grep Date_GenBank data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        LastUpdateDate=$(grep LastUpdateDate data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        FtpPath_RefSeq=$(grep FtpPath_RefSeq data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
        FtpPath_GenBank=$(grep FtpPath_GenBank data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | sed -e 's/^$/NA/' )
            
        total_length=$((cat data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | xtract -pattern DocumentSummary -element  BioSampleAccn -block Stat -if "@category" -equals total_length -element Stat) | awk -F'\t' '{print $2}' | head -1 | sed -e 's/^$/NA/')
        contig_count=$((cat data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | xtract -pattern DocumentSummary -element  BioSampleAccn -block Stat -if "@category" -equals contig_count -element Stat) | awk -F'\t' '{print $2}' | head -1 | sed -e 's/^$/NA/')
        contig_n50=$((cat data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | xtract -pattern DocumentSummary -element  BioSampleAccn -block Stat -if "@category" -equals contig_n50 -element Stat) | awk -F'\t' '{print $2}' | head -1 | sed -e 's/^$/NA/');

        #----------------------------------------------------------------------
        ## Assembly Biosample metadata
        esearch -db biosample -query "$AssemblyBiosampleAcc" | efetch -format native > data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab

        strain=$(grep "strain=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/strain=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        collection_date=$(grep "collection date=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/collection date=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        geographic_location=$(grep "geographic location=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/geographic location=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        isolation_source=$(grep "isolation source=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/isolation source=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        host=$(grep "host=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/host=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed -e 's/    //' )
        host_disease=$(grep "host disease=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/host disease=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        host_age=$(grep "host age=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/host age=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        collected_by=$(grep "collected by=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/collected by=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        serovar=$(grep "serovar=" data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed "s/\/serovar=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        biosample=$(sed -n -e 's/^.*BioSample: //p' data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab | sed 's/;//g' | awk '{print $1}' | head -1  | sed -e 's/^$/NA/' | sed 's/    //' )
        #----------------------------------------------------------------------

        echo -e "$AssemblyBiosampleAcc\t$strain\t$collection_date\t$geographic_location\t$isolation_source\t$host\t$host_disease\t$host_age\t$collected_by\t$serovar\t$AssemblyAccession\t$AssemblyName\t$AssemblyStatus\t$Isolate\t$Sub_value\t$WGS\t$total_length\t$contig_count\t$contig_n50\t$Taxid\t$Organism\t$SpeciesName\t$FtpPath_RefSeq\t$FtpPath_GenBank\t$SubmitterOrganization\t$SubmissionDate\t$AsmReleaseDate_GenBank\t$Date_GenBank\t$LastUpdateDate" > data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.2.tmp
    done
    sed -i 's/ /_/g' data/ncbi/metadata/tmp/*.Assembly-BioSampleAccn.2.tmp
    cat data/ncbi/metadata/tmp/*.Assembly-BioSampleAccn.2.tmp > data/ncbi/metadata/update.assemblies.metadata.tab
    sed -i -e '1i BioSampleAccn\tstrain\tcollection_date\tgeographic_location\tisolation_source\thost\thost_disease\thost_age\tcollected_by\tserovar\tAssemblyAccession\tAssemblyName\tAssemblyStatus\tIsolate\tstrain(Sub_value)\tWGS\ttotal_length\tcontig_count\tcontig_n50\tTaxid\tOrganism\tSpeciesName\tFtpPath_RefSeq\tFtpPath_GenBank\tSubmitterOrganization\tSubmissionDate\tAsmReleaseDate_GenBank\tDate_GenBank\tLastUpdateDate' data/ncbi/metadata/update.assemblies.metadata.tab
###############################################################################
## Download the files and arrange in 04_assembly
    if [ $download_ans == "y" ] ; then
        (mkdir results) > /dev/null 2>&1
        (mkdir results/04_assembly) > /dev/null 2>&1
        (mkdir results/04_assembly/all_fasta) > /dev/null 2>&1
        (mkdir data/ncbi/refseq) > /dev/null 2>&1
        for BioSample in $(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt ) ; do
            if [ ! -f results/04_assembly/all_fasta/$BioSample.fasta ] ; then
                echo "downloading $BioSample"

                (mkdir data/ncbi/refseq/$BioSample) > /dev/null 2>&1
                (mkdir results/04_assembly/$BioSample) > /dev/null 2>&1
                #esearch -db nucleotide -query $BioSample | efetch -format fasta > data/ncbi/refseq/$BioSample/$BioSample.fna
                #my_accession=$(esearch -db nucleotide -query $BioSample | efetch -format docsum | grep 'AssemblyAcc' | head -1 | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | head -1 )
                #esearch -db nucleotide -query $my_accession | efetch -format fasta > data/ncbi/refseq/$BioSample/$BioSample.fna
                ## this code was used to download samples that only appear in Biosample but not in Sra and assembly

                link1=$(esearch -db biosample -query "$BioSample" | elink -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank)
                link2=$(echo $link1 | awk -F'/' '{print $NF}')
                link=$(echo $link1 "/" $link2 "_genomic.fna.gz" | sed 's/ //g' )
                #echo $link
                (wget $link -P data/ncbi/refseq/$BioSample/) > /dev/null 2>&1
                gunzip -c data/ncbi/refseq/$BioSample/*.gz > data/ncbi/refseq/$BioSample/$BioSample.fna

                echo "arranging files in results/04_assembly"
                cp data/ncbi/refseq/$BioSample/$BioSample.fna results/04_assembly/$BioSample/$BioSample.fna
                cp results/04_assembly/$BioSample/$BioSample.fna results/04_assembly/$BioSample/contigs.fasta
                cp results/04_assembly/$BioSample/contigs.fasta results/04_assembly/$BioSample/$BioSample.fasta
                cp results/04_assembly/$BioSample/$BioSample.fasta results/04_assembly/$BioSample/$BioSample.joined.fasta
                sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/$BioSample/$BioSample.joined.fasta
                sed -i "1i "'>'$BioSample"" results/04_assembly/$BioSample/$BioSample.joined.fasta
                cp results/04_assembly/$BioSample/$BioSample.joined.fasta results/04_assembly/all_fasta/$BioSample.fasta
                sed -i '/^[[:space:]]*$/d' results/04_assembly/all_fasta/$BioSample.fasta
                fasta_formatter -i results/04_assembly/all_fasta/$BioSample.fasta -w 60 -o results/04_assembly/all_fasta/$BioSample.fasta.tmp
                mv results/04_assembly/all_fasta/$BioSample.fasta.tmp results/04_assembly/all_fasta/$BioSample.fasta
                else
                echo "$BioSample already downloaded"
            fi
        done
    fi
###############################################################################
## for type strains
    if [ "$type_ans" == "y" ] ; then
        (mkdir results/04_assembly/types_strain ) > /dev/null 2>&1
        for BioSample in $(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt ) ; do 
            name_tmp=$(grep $BioSample data/ncbi/tmp/all.accession.biosample.tab | awk '{print $4}' | awk -F'_' '{print $1, $2}' | sed 's/ /_/g' | wc -c )
            if [ $name_tmp -gt "2" ] ; then 
                name_1=$(grep $BioSample data/ncbi/tmp/all.accession.biosample.tab | awk '{print $4}' | awk -F'_' '{print $1, $2}' | sed 's/ /_/g' )
                name_2=$(grep $BioSample data/ncbi/tmp/all.accession.biosample.tab | awk '{print $5}' | sed 's/ /_/g' | sed 's/$/_T/g' )
                name=$(echo $name_1"_"$name_2 )
                cp results/04_assembly/all_fasta/$BioSample.fasta results/04_assembly/types_strain/$name.fasta
                else
                name=$(echo $BioSample)
                cp results/04_assembly/all_fasta/$BioSample.fasta results/04_assembly/types_strain/$name.fasta
            fi
            sed -i 's/>.*/>'$name'/' results/04_assembly/types_strain/$name.fasta
        done
    fi
###############################################################################
