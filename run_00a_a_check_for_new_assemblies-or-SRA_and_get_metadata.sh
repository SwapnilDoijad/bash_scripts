#!/bin/bash
###############################################################################
## DO NOT DELETE metadata/tmp. VERY VERY IMP.
###############################################################################
echo "Hi Swapnil, try to modify the script so you can have data for Pacbio and Nanopore \
the final output should be (after checking duplicates) \
data/ncbi/metadata/update.assemblies.metadata.tab "
###############################################################################
# GENOMIC PAIRED ILLUMINA
# GENOMIC SINGLE ILLUMINA
# GENOMIC SINGLE LS454
# GENOMIC SINGLE OXFORD_NANOPORE
# GENOMIC SINGLE PACBIO_SMRT
# METAGENOMIC PAIRED ILLUMINA
# METAGENOMIC SINGLE ION_TORRENT
# OTHER PAIRED ILLUMINA
# OTHER SINGLE PACBIO_SMRT
# TRANSCRIPTOMIC PAIRED ILLUMINA
# TRANSCRIPTOMIC SINGLE ILLUMINA
# TRANSCRIPTOMIC SINGLE PACBIO_SMRT
###############################################################################
## Initial input and directory creation 
    PATH=$PATH:$HOME/tools/edirect
    T=`date '+%d%m%Y_%H%M%S'`

    (mkdir data) > /dev/null 2>&1
    (mkdir data/ncbi) > /dev/null 2>&1
    (mkdir data/ncbi/tmp) > /dev/null 2>&1
    (mkdir data/ncbi/metadata) > /dev/null 2>&1
    (mkdir data/ncbi/metadata/tmp) > /dev/null 2>&1
    #------------------------------------------------------------------------------
    ## genomic Prf details
    ## presently only going designed for illumina but easily extended to other  
    LS=$(echo "GENOMIC") #LibrarySource
    LL=$(echo "SINGLE") #LibraryLayout
    Plf=$(echo "PACBIO_SMRT") #Platform
    #------------------------------------------------------------------------------
    echo “Update: What is the name of the genus species?”
    read name
    echo "Wish to check assembly-SRA newly? answer y or PRESS ENTER if already checked" 
    read assemblySRA_check
    echo "Do you want to get the metadata for new assemblies and SRA? answer a for assemblies, s for SRA and as for both" 
    read metadata
###############################################################################
## Checking new Assembly-BioSampleAccn
    if [ "$assemblySRA_check" == "y" ]; then
    #echo "Checking new Assembly-BioSampleAccn for $name"
    ## find all the BioSamples
    (esearch -db assembly -query "$name" | efetch -format docsum | xtract -pattern DocumentSummary -element BioSampleAccn) > data/ncbi/tmp/"$name".BioSampleAccn."$T".tmp
    total_assemblies=$(cat data/ncbi/tmp/"$name".BioSampleAccn."$T".tmp | awk '{print $1}' | sort -u | wc -l)
    ## check if these not present previously
    (ls data/ncbi/metadata/tmp/*.Assembly-BioSampleAccn.tmp | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' | sort > data/ncbi/metadata/tmp/old_BioSampleAccn.txt) > /dev/null 2>&1
    (cat data/ncbi/tmp/"$name".BioSampleAccn."$T".tmp | sort -u > data/ncbi/metadata/tmp/new_BioSampleAccn.txt)> /dev/null 2>&1
    comm -13 data/ncbi/metadata/tmp/old_BioSampleAccn.txt data/ncbi/metadata/tmp/new_BioSampleAccn.txt > data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt 
    cp data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt data/ncbi/tmp/list."$name".Assembly-BioSampleAccn.to-update."$T".txt 
    rm data/ncbi/tmp/"$name".BioSampleAccn."$T".tmp
    rm data/ncbi/metadata/tmp/old_BioSampleAccn.txt
    rm data/ncbi/metadata/tmp/new_BioSampleAccn.txt
    total_assemblies_new=$(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt | awk '{print $1}' | sort -u | wc -l)
    ## get all-BioSampleAcc (I am not sure if this is covered already)
    cat data/ncbi/metadata/tmp/old_BioSampleAccn.txt data/ncbi/metadata/tmp/new_BioSampleAccn.txt | sort | uniq > data/ncbi/list."$name".Assembly-BioSampleAccn.all-BioSampleAccn.txt
###############################################################################
## Checking new SRA-BioSample 
    #echo "Checking new SRA-BioSample for $name"
    ## find all the BioSamples
    (esearch -db sra -query "$name" | esummary -format runinfo -mode xml | xtract -pattern Row -element BioSample,LibrarySource,LibraryLayout,Platform) > data/ncbi/tmp/"$name".SRA-BioSample."$T".all.tmp
    ## stat
    echo "-------------------------------------------------------------------------"
    awk '{first = $1; $1 = ""; print $0 }' data/ncbi/tmp/"$name".SRA-BioSample."$T".all.tmp | sort | uniq -c
    echo "-------------------------------------------------------------------------"
    awk '{first = $1; $1 = ""; print $0 }' data/ncbi/tmp/"$name".SRA-BioSample."$T".all.tmp | sort | uniq -c > data/ncbi/"$name".SRA-BioSample."$T".all.tab
    ## short-list according to LS LL Prf
    awk -F'\t' '($2 == "'$LS'" && $3 == "'$LL'" && $4 == "'$Plf'" ) {print $1}' data/ncbi/tmp/"$name".SRA-BioSample."$T".all.tmp > data/ncbi/tmp/"$name".SRA-BioSample.$LS.$LL.$Plf."$T".tmp
    total_SRA=$(cat data/ncbi/tmp/"$name".SRA-BioSample.$LS.$LL.$Plf."$T".tmp | awk '{print $1}' | sort -u | wc -l)
    ## check if these not present previously
    (ls data/ncbi/metadata/tmp/*_SRA.metadata.$LS.$LL.$Plf.tab | awk -F'/' '{print $5}' | awk -F'_' '{print $1}' | sort > data/ncbi/metadata/tmp/old_BioSample.$LS.$LL.$Plf.tmp) > /dev/null 2>&1
    ## remove samples those are already identified as duplicates in assembly in last run
    ## create file ..as_well.txt else for first time, error will appear on screen
    if [ ! -f data/ncbi/tmp/list."$name".SRA-BioSample.duplicate.in_assembly_as_well.$LS.$LL.$Plf.txt ] ; then touch data/ncbi/tmp/list."$name".SRA-BioSample.duplicate.in_assembly_as_well.$LS.$LL.$Plf.txt ; fi
    (cat data/ncbi/tmp/list."$name".SRA-BioSample.duplicate.in_assembly_as_well.$LS.$LL.$Plf.txt data/ncbi/metadata/tmp/old_BioSample.$LS.$LL.$Plf.tmp | sort > data/ncbi/metadata/tmp/old_BioSample.$LS.$LL.$Plf.txt) > /dev/null 2>&1
    (cat data/ncbi/tmp/"$name".SRA-BioSample.$LS.$LL.$Plf."$T".tmp | sort -u > data/ncbi/metadata/tmp/new_BioSample.$LS.$LL.$Plf.txt) > /dev/null 2>&1
    comm -13 data/ncbi/metadata/tmp/old_BioSample.$LS.$LL.$Plf.txt data/ncbi/metadata/tmp/new_BioSample.$LS.$LL.$Plf.txt > data/ncbi/list."$name".SRA-BioSample.to-update.$LS.$LL.$Plf.txt
    cp data/ncbi/list."$name".SRA-BioSample.to-update.$LS.$LL.$Plf.txt data/ncbi/tmp/list."$name".SRA-BioSample.to-update."$T".$LS.$LL.$Plf.txt
    rm data/ncbi/tmp/"$name".SRA-BioSample.$LS.$LL.$Plf."$T".tmp
    rm data/ncbi/metadata/tmp/old_BioSample.$LS.$LL.$Plf.txt
    rm data/ncbi/metadata/tmp/new_BioSample.$LS.$LL.$Plf.txt
    rm data/ncbi/metadata/tmp/old_BioSample.$LS.$LL.$Plf.tmp
    total_SRA_new=$(cat data/ncbi/list."$name".SRA-BioSample.to-update.$LS.$LL.$Plf.txt | awk '{print $1}' | sort -u | wc -l )
###############################################################################
## remove duplicate Assembly-BioSampleAccn from SRA-BioSample
    #echo "checking duplicates in Assembly-BioSampleAccn from SRA-BioSample"
        for F1 in $(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt); do
            V3=$(grep $F1 data/ncbi/list."$name".SRA-BioSample.to-update.$LS.$LL.$Plf.txt)
            if [ ! -z $V3 ]; then ## if SRA accession not in Assembly-BioSampleAccn then
                sed -i "/$F1/d" data/ncbi/list."$name".SRA-BioSample.to-update.$LS.$LL.$Plf.txt ## delete from SRA to-update list
                echo $F1 >> data/ncbi/tmp/list."$name".SRA-BioSample.duplicate.in_assembly_as_well.$LS.$LL.$Plf.txt ## write to duplicate list
            fi
        done
    ##-------------------------------------------------------------------------
    ## If the biosample is absent in the assembly, check short reads and long reads, both, if present for a biosample
        cat data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.PAIRED.ILLUMINA.txt | sort > data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.PAIRED.ILLUMINA.sorted.txt
        cat data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.SINGLE.OXFORD_NANOPORE.txt data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.SINGLE.PACBIO_SMRT.txt | sort > \
        data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.SINGLE.long_read.sorted.txt
        comm -12 data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.PAIRED.ILLUMINA.sorted.txt data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.SINGLE.long_read.sorted.txt > \
        data/ncbi/list."$name".SRA-BioSample.short_long_read_both_available.txt
        rm data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.PAIRED.ILLUMINA.sorted.txt
        rm data/ncbi/list."$name".SRA-BioSample.to-update.GENOMIC.SINGLE.long_read.sorted.txt
    ##-------------------------------------------------------------------------

    new_assemblies=$(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt | awk '{print $1}' | sort -u | wc -l )
    new_SRA=$(cat data/ncbi/list."$name".SRA-BioSample.to-update.$LS.$LL.$Plf.txt | awk '{print $1}' | sort -u | wc -l )
    duplicates=$(cat data/ncbi/tmp/list."$name".SRA-BioSample.duplicate.in_assembly_as_well.$LS.$LL.$Plf.txt | awk '{print $1}' | sort -u | wc -l )
    #echo "-------------------------------------------------------------------------"
    echo "total assemblies available:            $total_assemblies"
    echo "total assemblies available NEW:        $total_assemblies_new"
    echo "total SRA available:                   $total_SRA"
    echo "total SRA available NEW:               $total_SRA_new"
    echo "SRA (already asmbly), not considering: $duplicates"
    echo "new Assemblies now to process:         $new_assemblies"
    echo "new SRA now to process:                $new_SRA ($LS $LL $Plf)"
    echo "-------------------------------------------------------------------------"
    fi
###############################################################################
## metadata
    if [ "$metadata" == "as" ] || [ "$metadata" == "a" ] ; then
    #------------------------------------------------------------------------------
    ## assembly metadata
        for AssemblyBiosampleAcc in $(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt); do
            if [ -s data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp ] && [ -s data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.2.tmp ]; then
            echo "metadata for $AssemblyBiosampleAcc assembly already extracted"
            else
            echo "extracting Assembly-BioSampleAccn metadata for $AssemblyBiosampleAcc"
            esearch -db assembly -query "$AssemblyBiosampleAcc" | efetch -format docsum > data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp
            ## check if the data is extracted (sometimes NCBI server fails)
            if [ -s data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp ]; then
                #echo "collecting data" 
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
                contig_n50=$((cat data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.tmp | xtract -pattern DocumentSummary -element  BioSampleAccn -block Stat -if "@category" -equals contig_n50 -element Stat) | awk -F'\t' '{print $2}' | head -1 | sed -e 's/^$/NA/')
                #----------------------------------------------------------------------
                ## Assembly Biosample metadata
                esearch -db biosample -query "$AssemblyBiosampleAcc" | efetch -format native > data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab
                if [ -s data/ncbi/metadata/tmp/$AssemblyBiosampleAcc.tmp.tab ] ; then 
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
                    else
                    echo "esearch for biosample failing"
                fi
                #----------------------------------------------------------------------
                echo -e "$AssemblyBiosampleAcc\t$strain\t$collection_date\t$geographic_location\t$isolation_source\t$host\t$host_disease\t$host_age\t$collected_by\t$serovar\t$AssemblyAccession\t$AssemblyName\t$AssemblyStatus\t$Isolate\t$Sub_value\t$WGS\t$total_length\t$contig_count\t$contig_n50\t$Taxid\t$Organism\t$SpeciesName\t$FtpPath_RefSeq\t$FtpPath_GenBank\t$SubmitterOrganization\t$SubmissionDate\t$AsmReleaseDate_GenBank\t$Date_GenBank\t$LastUpdateDate" > data/ncbi/metadata/tmp/"$AssemblyBiosampleAcc".Assembly-BioSampleAccn.2.tmp
            fi
            fi
        done

        sed -i 's/ /_/g' data/ncbi/metadata/tmp/*.Assembly-BioSampleAccn.2.tmp
        cat data/ncbi/metadata/tmp/*.Assembly-BioSampleAccn.2.tmp > data/ncbi/metadata/update.assemblies.metadata.tab
        sed -i -e '1i BioSampleAccn\tstrain\tcollection_date\tgeographic_location\tisolation_source\thost\thost_disease\thost_age\tcollected_by\tserovar\tAssemblyAccession\tAssemblyName\tAssemblyStatus\tIsolate\tstrain(Sub_value)\tWGS\ttotal_length\tcontig_count\tcontig_n50\tTaxid\tOrganism\tSpeciesName\tFtpPath_RefSeq\tFtpPath_GenBank\tSubmitterOrganization\tSubmissionDate\tAsmReleaseDate_GenBank\tDate_GenBank\tLastUpdateDate' data/ncbi/metadata/update.assemblies.metadata.tab
        fi
    ###############################################################################
    ## SRA-metadata 
    if [ "$metadata" == "as" ] || [ "$metadata" == "s" ]; then
        for F1 in $(cat data/ncbi/list."$name".SRA-BioSample.to-update.$LS.$LL.$Plf.txt); do
        if [ -s data/ncbi/metadata/tmp/"$F1"_SRA.metadata.$LS.$LL.$Plf.tab ]; then
            echo "metadata for $F1 $LS.$LL.$Plf already extracted"
        else
            (esearch -db sra -query "$F1" | esummary -format runinfo -mode xml | xtract -pattern Row -element BioSample,Run,bases,avgLength,size_MB,Experiment,LibrarySource,LibraryLayout,Platform,Model,SRAStudy,BioProject,ProjectID,Sample,SampleType,TaxID,ScientificName,SampleName,Submission,ReleaseDate,download_path) | head -1 > data/ncbi/metadata/tmp/"$F1"_SRA.metadata.$LS.$LL.$Plf.tab
            #------------------------------------------------------------------------------
            ## SRA--biosample-metadata: input SRA, links to biosample and get all possible meadata (then needs to define what needs to be extracted for e.g strain, host, date etc.)

            echo "extracting SRA-BioSample metadata for $F1 for $LS $LL $Plf"

            awk -F'\t' '($7 == "'$LS'" && $8 == "'$LL'" && $9 == "'$Plf'" ) {print $14}'  data/ncbi/metadata/tmp/"$F1"_SRA.metadata.$LS.$LL.$Plf.tab > data/ncbi/tmp/$F1.list.biosample.accessions.$LS.$LL.$Plf.txt
        ##
            for F3 in $(cat data/ncbi/tmp/$F1.list.biosample.accessions.$LS.$LL.$Plf.txt); do

                esearch -db biosample -query "$F3" | efetch -format native > data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab

                strain=$(grep "strain=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/strain=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                collection_date=$(grep "collection date=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/collection date=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                geographic_location=$(grep "geographic location=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/geographic location=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                isolation_source=$(grep "isolation source=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/isolation source=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                host=$(grep "host=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/host=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                host_disease=$(grep "host disease=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/host disease=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                host_age=$(grep "host age=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/host age=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                collected_by=$(grep "collected by=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/collected by=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                serovar=$(grep "serovar=" data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed "s/\/serovar=//g" | sed "s/\"//g" | head -1  | sed -e 's/^$/NA/')
                biosample=$(sed -n -e 's/^.*BioSample: //p' data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed 's/;//g' | awk '{print $1}' | head -1  | sed -e 's/^$/NA/')
                SRA=$(sed -n -e 's/^.*SRA: //p' data/ncbi/metadata/tmp/$F1.$F3.tmp.$LS.$LL.$Plf.tab | sed 's/;//g' | awk '{print $1}' | head -1  | sed -e 's/^$/NA/')

                echo -e "$biosample\t$strain\t$collection_date\t$geographic_location\t$isolation_source\t$host\t$host_disease\t$host_age\t$collected_by\t$serovar\t$SRA" > data/ncbi/metadata/tmp/"$F1"_SRA.metadata.metadata.$LS.$LL.$Plf.tab
                sed -i "s/    //g" data/ncbi/metadata/tmp/"$F1"_SRA.metadata.metadata.$LS.$LL.$Plf.tab

                if [ ! -s data/ncbi/metadata/tmp/"$F1"_SRA.metadata.metadata.$LS.$LL.$Plf.tab ] ; then
                echo "data collection for $F1 $LS.$LL.$Plf failed"
                rm -rf data/ncbi/metadata/tmp/"$F1"_SRA.metadata.metadata.$LS.$LL.$Plf.tab
                fi
            done
            rm data/ncbi/tmp/*.list.biosample.accessions.$LS.$LL.$Plf.txt
        fi
        done
        #------------------------------------------------------------------------------
        cat data/ncbi/metadata/tmp/*SRA.metadata.$LS.$LL.$Plf.tab > data/ncbi/metadata/tmp/update.SRA.metadata.$LS.$LL.$Plf.tab1
        sed -i -e '1i	BioSample\tRun\tbases\tavgLength\tsize_MB\tExperiment\tLibrarySource\tLibraryLayout\tPlatform\tModel\tSRAStudy\tBioProject\tProjectID\tSample\tSampleType\tTaxID\tScientificName\tSampleName\tSubmission\tReleaseDate\tdownload_path' data/ncbi/metadata/tmp/update.SRA.metadata.$LS.$LL.$Plf.tab1
        cat data/ncbi/metadata/tmp/*SRA.metadata.metadata.$LS.$LL.$Plf.tab > data/ncbi/metadata/tmp/update.SRA.metadata.metadata.$LS.$LL.$Plf.tab1 ;
        sed -i -e '1i	BioSample\tstrain\tcollection_date\tgeographic_location\tisolation_source\thost\thost_disease\thost_age\tcollected_by\tserovar\tSRA' data/ncbi/metadata/tmp/update.SRA.metadata.metadata.$LS.$LL.$Plf.tab1
        #------------------------------------------------------------------------------
        paste data/ncbi/metadata/tmp/update.SRA.metadata.metadata.$LS.$LL.$Plf.tab1 data/ncbi/metadata/tmp/update.SRA.metadata.$LS.$LL.$Plf.tab1 > data/ncbi/metadata/update.SRA.metadata.$LS.$LL.$Plf.tab
    fi
###############################################################################
