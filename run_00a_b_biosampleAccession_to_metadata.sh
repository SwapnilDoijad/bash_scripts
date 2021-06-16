(mkdir data) > /dev/null 2>&1
(mkdir data/ncbi) > /dev/null 2>&1
(mkdir data/ncbi/tmp) > /dev/null 2>&1
(mkdir data/ncbi/metadata) > /dev/null 2>&1
(mkdir data/ncbi/metadata/tmp) > /dev/null 2>&1
## assembly metadata
    #for AssemblyBiosampleAcc in $(cat data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt); do
    echo "provide the name of the biosamples list (for e.g. list.biosamples.txt )"
	read list
    echo "name of the genus and species (for e.g. Escherichia coli)"
    read name
    for AssemblyBiosampleAcc in $(cat $list); do
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
    awk 'NR>1 {print $1}' data/ncbi/metadata/update.assemblies.metadata.tab > data/ncbi/list."$name".Assembly-BioSampleAccn.to-update.txt
