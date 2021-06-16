## Initial input and directory creation 
    PATH=$PATH:$HOME/tools/edirect
    T=`date '+%d%m%Y_%H%M%S'`

    (mkdir data) > /dev/null 2>&1
    (mkdir data/ncbi) > /dev/null 2>&1
    (mkdir data/ncbi/tmp) > /dev/null 2>&1
    (mkdir data/ncbi/metadata_Abr) > /dev/null 2>&1
    (mkdir data/ncbi/metadata_Abr/raw_files) > /dev/null 2>&1
    (mkdir data/ncbi/metadata_Abr/raw_files/data_was_not_available) > /dev/null 2>&1
    (mkdir data/ncbi/metadata_Abr/tmp) > /dev/null 2>&1

for biosample in $(cat list.all.biosamples.txt); do
    ## check if phepotype if already extracted
    if [ -f data/ncbi/metadata_Abr/raw_files/$biosample.Abr_phenotype.tab ] ; then #|| [ ! -s data/ncbi/metadata_Abr/raw_files/data_was_not_available/$biosample.Abr_phenotype.tab ]; then 
        echo "Abr phenotype for $biosample is already present"
        else
        echo "finding Abr phenotype for $biosample"
        esearch -db biosample -query "$biosample" | efetch -format docsum | xtract -pattern DocumentSummary -group Comment -block Row -sep "|" -element Cell | sed 's/\t/\n/g' | sed 's/|/\t/g' > data/ncbi/metadata_Abr/raw_files/$biosample.Abr_phenotype.tab
        ## if file is empty move file to data_was_not_available folder. This will ease re-running the program. 
        ## if you want to check if data if newly added then remove second condition above
        if [ ! -s data/ncbi/metadata_Abr/raw_files/$biosample.Abr_phenotype.tab ]  ; then
            mv data/ncbi/metadata_Abr/raw_files/$biosample.Abr_phenotype.tab data/ncbi/metadata_Abr/raw_files/$biosample.Abr_phenotype.tab
        fi
    fi
done
sed -i 's/ /_/g' data/ncbi/metadata_Abr/raw_files/*.Abr_phenotype.tab

###############################################################################
## Antibiotics and number of isolates for which data available
echo "calculating Antibiotics and number of isolates for which data available"
(mkdir data/ncbi/metadata_Abr/tmp) > /dev/null 2>&1

no_of_isolates=$(ls data/ncbi/metadata_Abr/raw_files/*.tab | wc -l)
awk '{print $1}' data/ncbi/metadata_Abr/raw_files/*.tab | sort -u > data/ncbi/metadata_Abr/list.ab.txt
cat data/ncbi/metadata_Abr/raw_files/*.tab | awk -F'\t' '{print $2}' | sort -u > data/ncbi/metadata_Abr/tmp/pattern.txt

(rm data/ncbi/metadata_Abr/tmp/ab_class.tab) > /dev/null 2>&1
for Antb in $(cat data/ncbi/metadata_Abr/list.ab.txt); do
    data_available_for=$(grep $Antb data/ncbi/metadata_Abr/raw_files/*.tab | awk '{print $1}' | uniq | awk -F':' '$2=="'$Antb'" {print $2}' | wc -l)
    (rm data/ncbi/metadata_Abr/tmp/$Antb.tmp) > /dev/null 2>&1
    for pattern in $(cat data/ncbi/metadata_Abr/tmp/pattern.txt); do
        total=$(grep $Antb data/ncbi/metadata_Abr/raw_files/*.tab | awk -F'\t' '$2 == "'$pattern'" {print $0}' | awk '{print $1}' | uniq | awk -F':' '$2=="'$Antb'" {print $2}' | wc -l )
        percnt=$( echo 'scale=1;100*'$total'/'$data_available_for'' | bc ) 
        echo $pattern $total $percnt >> data/ncbi/metadata_Abr/tmp/$Antb.tmp
    done
    pattern_2=$(cat data/ncbi/metadata_Abr/tmp/$Antb.tmp | awk '{print $2, $3}' | tr "\n" " " )
    pattern_2_sum=$(awk '{sum += $2} END {print sum}' data/ncbi/metadata_Abr/tmp/$Antb.tmp)
    variables=$(( $data_available_for - $pattern_2_sum ))
    echo $Antb $no_of_isolates $data_available_for $pattern_2 $variables >> data/ncbi/metadata_Abr/antibiotics_phenotypes.$T.csv
    (rm data/ncbi/metadata_Abr/tmp/*.tmp) > /dev/null 2>&1
    
    #--------------
    ab_class=$(grep -i "$Antb" /media/network/gene_databases/ab_class.tab )
    if [ -z "$ab_class" ] ; then 
        ab_class_2=$(echo "unknown")
        else
        ab_class_2=$(echo "$ab_class" | awk -F'\t' '{$1=""; print $0}' | uniq | sed 's/ /_/g' | sed -e "s/^_//" | tr "\n" "_" | tr " " "_" | sed 's/_$//' )
    fi
    echo $Antb $ab_class_2 >> data/ncbi/metadata_Abr/tmp/ab_class.tab
done
sed -i '1 i\ab_class ab_class' data/ncbi/metadata_Abr/tmp/ab_class.tab

sed -i 's/$/ \%/g' data/ncbi/metadata_Abr/tmp/pattern.txt
pattern_3=$(cat data/ncbi/metadata_Abr/tmp/pattern.txt | tr "\n" " " )
echo Antibiotic no_of_isolates data_available_for $pattern_3 variables | cat - data/ncbi/metadata_Abr/antibiotics_phenotypes.$T.csv > temp && mv temp data/ncbi/metadata_Abr/antibiotics_phenotypes.$T.csv
paste data/ncbi/metadata_Abr/tmp/ab_class.tab data/ncbi/metadata_Abr/antibiotics_phenotypes.$T.csv > temp && mv temp data/ncbi/metadata_Abr/antibiotics_phenotypes.$T.csv 

###############################################################################
## list out all/susceptible/resistant isolates
echo "creating list for susceptible_resistant"
(mkdir data/ncbi/metadata_Abr/susceptible_resistant)> /dev/null 2>&1
(mkdir data/ncbi/metadata_Abr/susceptible_resistant/tmp)> /dev/null 2>&1
(rm data/ncbi/metadata_Abr/susceptible_resistant/traits.dbgwas.$Antb.txt)> /dev/null 2>&1
for Antb in $(cat data/ncbi/metadata_Abr/list.ab.txt); do
## get list ALL
grep -i "$Antb" data/ncbi/metadata_Abr/raw_files/*.tab | sed 's/data\/ncbi\/metadata_Abr\/raw_files\///g' | sed -e 's/\.Abr_phenotype\.tab\:/\t/g' | sed 's/ /_/g' > data/ncbi/metadata_Abr/susceptible_resistant/tmp/$Antb.tab
## get list resistant
grep -i "$Antb" data/ncbi/metadata_Abr/raw_files/*.tab | awk -F'\t' '$2 == "resistant" {print $0}' | awk '{print $1, $2}' | uniq | awk -F':' '$2=="'$Antb' resistant" {print $0}' | sed 's/data\/ncbi\/metadata_Abr\/raw_files\///g' | sed -e 's/\.Abr_phenotype\.tab\:/\t/g' | sed 's/ /_/g' > data/ncbi/metadata_Abr/susceptible_resistant/tmp/$Antb.resistant.tab
## get list susceptible
grep -i "$Antb" data/ncbi/metadata_Abr/raw_files/*.tab | awk -F'\t' '$2 == "susceptible" {print $0}' | awk '{print $1, $2}' | uniq | awk -F':' '$2=="'$Antb' susceptible" {print $0}' | sed 's/data\/ncbi\/metadata_Abr\/raw_files\///g' | sed -e 's/\.Abr_phenotype\.tab\:/\t/g' | sed 's/ /_/g' > data/ncbi/metadata_Abr/susceptible_resistant/tmp/$Antb.susceptible.tab
awk '{print $1",1"}' data/ncbi/metadata_Abr/susceptible_resistant/tmp/$Antb.resistant.tab > data/ncbi/metadata_Abr/susceptible_resistant/traits.dbgwas.$Antb.txt
awk '{print $1",0"}' data/ncbi/metadata_Abr/susceptible_resistant/tmp/$Antb.susceptible.tab >> data/ncbi/metadata_Abr/susceptible_resistant/traits.dbgwas.$Antb.txt

## remove duplicates (isolates with resistant as well as susceptible phenotype)
awk -F',' '{print $1}' data/ncbi/metadata_Abr/susceptible_resistant/traits.dbgwas.$Antb.txt | sort | uniq -d > data/ncbi/metadata_Abr/susceptible_resistant/tmp/susceptible_resistant.txt
    for dupl in $(cat data/ncbi/metadata_Abr/susceptible_resistant/tmp/susceptible_resistant.txt) ; do
    sed -i '/'$dupl'/d' data/ncbi/metadata_Abr/susceptible_resistant/traits.dbgwas.$Antb.txt
    done
done

###############################################################################
## summarise all the metadata downloaded
(rm data/ncbi/metadata.stat.csv)> /dev/null 2>&1
for Antb in $(cat data/ncbi/metadata_Abr/list.ab.txt); do
total=$(cat data/ncbi/metadata_Abr/susceptible_resistant/traits.dbgwas.$Antb.txt | wc -l)
R=$(cat data/ncbi/metadata_Abr/susceptible_resistant/traits.dbgwas.$Antb.txt | awk -F',' '{print $2}' | grep -c '1')
S=$(cat data/ncbi/metadata_Abr/susceptible_resistant/traits.dbgwas.$Antb.txt | awk -F',' '{print $2}' | grep -c '0')
echo $Antb $total $R $S >> data/ncbi/metadata.stat.csv
done
sed  -i '1i Antb total R S' data/ncbi/metadata.stat.csv
###############################################################################
