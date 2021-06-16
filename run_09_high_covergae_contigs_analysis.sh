#!/bin/bash
###############################################################################
# high coverage contigs

###############################################################################

echo "started.... step-09 high coverage contigs analysis ----------------------"

###############################################################################
echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

#echo "want to suffix files? type "y" else press enter to continue"
#read answer
#if [ "$answer" = "y" ]; then
#s="_$l"
#fi

###############################################################################

if [ ! -d results/09_high_coverage_contigs ]; then
mkdir results/09_high_coverage_contigs
fi

if [ ! -d results/09_high_coverage_contigs/tmp ]; then
mkdir results/09_high_coverage_contigs/tmp
fi

###############################################################################
## preliminary analysis for high coverage contig analysis

for F1 in $(cat $list); do

echo "running.... preliminary analysis for high coverage contigs for $F1"

if [ ! -d results/09_high_coverage_contigs/tmp/$F1 ]; then
mkdir results/09_high_coverage_contigs/tmp/$F1
fi

grep -F ">" results/04_assembly/raw_files/$F1/$F1.fasta | sed -e 's/_/ /g' | sort -nrk 6 | awk '$6>=10.0 && $4>=500 {print $0}' | sed -s 's/ /_/g' | sed -e 's/>//g' | awk -F'_' '{$6 = sprintf("%04d", $6); print}' | sort -nrk 6 | sed -e 's/ /_/g' > results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv

AverAssemContCov=$(grep ">" results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta | awk -F '_' '{ sum += $6; n++ } END { if (n > 0) print sum / n; }')

grep -F ">" results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta | sed -e 's/_/ /g' | sort -nrk 6 | awk '$6>=2*'$AverAssemContCov' && $4>=500 {print $0}' | sed -s 's/ /_/g' | sed -e 's/>//g' | awk -F'_' '{$6 = sprintf("%04d", $6); print}' | sort -nrk 6 | sed -e 's/ /_/g' > results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.2xHCv.csv

grep -F ">" results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta | sed -e 's/_/ /g' | sort -nrk 6 | awk '$6>=2*'$AverAssemContCov' && $4>=500 {print $0}' | sed -s 's/ /_/g' | sed -e 's/>//g' | awk -F'_' '{print $6}' | sort -nrk 6 | sed -e 's/ /_/g' > results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.2xHCv.csv.to_filter

perl  /home/swapnil/pipeline/tools/fastagrep.pl -f results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.2xHCv.csv.to_filter results/04_assembly/raw_files/$F1/$F1.fasta > results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta

cp results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta.joined.fasta 

sed -i 's/>.*/NNNNNNNNNN/g' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta.joined.fasta
sed -i "1i "'>'$F1"" results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta.joined.fasta

perl /home/swapnil/tools/prokka-1.13/bin/prokka --quiet --outdir results/09_high_coverage_contigs/tmp/$F1 --force --prefix $F1 --locustag $F1 --strain $F1  results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta.joined.fasta > /dev/null 2>&1

(ugene --task=/home/swapnil/pipeline/tools/gbk2csv-Ugene.uwl --in=results/09_high_coverage_contigs/tmp/$F1/$F1.gbk --out=results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv.tmp --format=csv) > /dev/null 2>&1
awk '$1 ~ /Group/ {print $0}' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv.tmp > results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv
awk '$1 ~ /CDS/ {print $0}' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv.tmp >> results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv

##------------------------------------------------------------------------------
## if the file doesnt contain translation column, add!

if grep -q translation results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv; then
cp results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.2.csv
else
awk '{print $0"SUFFIX"}' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv > results/09_high_coverage_contigs/tmp/$F1/tmp.csv
sed 's/SUFFIX/\,"translation"/g' results/09_high_coverage_contigs/tmp/$F1/tmp.csv > results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.2.csv
(rm tmp.csv) > /dev/null 2>&1
fi

##------------------------------------------------------------------------------

## separating each column
awk -F',' 'NR==1{for (i=1;i<=NF;i++) a[i]=$i; next} {for (i=1;i<=NF;i++) {print $i > a[i]".'$F1'.tmp"}}' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.2.csv ; rename 's/"//' *.$F1.tmp ; rename 's/"././' *.$F1.tmp
paste Start.$F1.tmp End.$F1.tmp Length.$F1.tmp locus_tag.$F1.tmp product.$F1.tmp translation.$F1.tmp >> results/09_high_coverage_contigs/all.HCv.csv

rm *.$F1.tmp
##------------------------------------------------------------------------------

echo "finished... preliminary analysis for high coverage contigs for $F1"
done
sed -i '1 i\Start,End,Length,locus_tag,product,translation' results/09_high_coverage_contigs/tmp/all.HCv.csv
ssconvert results/09_high_coverage_contigs/all.HCv.csv results/09_high_coverage_contigs/all.HCv.xlsx

echo "finished... writing the preliminary analysis results"

## FINISH: preliminary analysis for high coverage contig analysis
###############################################################################

###############################################################################
## separate high coverage contigs, annotate and write in csv

#------------------------------------------------------------------------------
#shorten contig name (so just NODE_XX, needed later for PROKKA)

for F1 in $(cat $list); do

echo "Started separatly writing high coverage contigs for $F1"

cp results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta.tmp
sed -i 's/_length.*//g' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta.tmp

#separately write each contigs to different file
while read line 
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=$F1.HCv.contig.${line#>}.fa
        echo $line > results/09_high_coverage_contigs/tmp/$F1/$outfile
    else
        echo $line >> results/09_high_coverage_contigs/tmp/$F1/$outfile
    fi
done < results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.contigs-filtered.fasta.tmp

cd results/09_high_coverage_contigs/tmp/$F1/
ls *.fa > list2.tmp
cd ..
cd ..
cd ..
cd ..
sed -i 's/\.fa//g' results/09_high_coverage_contigs/tmp/$F1/list2.tmp
cp results/09_high_coverage_contigs/tmp/$F1/list2.tmp results/09_high_coverage_contigs/tmp/$F1/list3.tmp
sed -i 's/'$F1'.HCv.contig.//g' results/09_high_coverage_contigs/tmp/$F1/list3.tmp

echo "Finished separatly writing high coverage contigs for $F1"

done

###############################################################################
## annotate high coverage separated contigs 
for F2 in $(cat $list); do
cd results/09_high_coverage_contigs/tmp/$F2

	for F3 in $(cat list2.tmp); do

	echo "Started anotating high coverage separated contigs for $F3"

	( /home/swapnil/tools/prokka-1.13/bin/prokka --quiet --outdir "$F3"_tmp --force --prefix $F3 --locustag $F3 --strain $F3 --rnammer $F3.fa) > /dev/null 2>&1

	grep "product=" "$F3"_tmp/$F3.gbk > "$F3"_tmp/$F3.gbk.txt 
	sed -i 's/ /_/g' "$F3"_tmp/$F3.gbk.txt
	sed -i 's/\/product="//g' "$F3"_tmp/$F3.gbk.txt
	sed -i 's/_____________________//g' "$F3"_tmp/$F3.gbk.txt
	sed -i 's/\"//g' "$F3"_tmp/$F3.gbk.txt


	awk '
	{ 
	    for (i=1; i<=NF; i++)  {
	        a[NR,i] = $i
	    }
	}
	NF>p { p = NF }
	END {    
	    for(j=1; j<=p; j++) {
	       str=a[1,j]
	       for(i=2; i<=NR; i++){
	           str=str" "a[i,j];
	       }
	       print str
	   }
	}' "$F3"_tmp/$F3.gbk.txt > "$F3"_tmp/$F3.gbk.txt.tmp ; sed -i 's/ /==/g' "$F3"_tmp/$F3.gbk.txt.tmp
	done

cd ..
cd ..
cd ..
cd ..

echo "Finished anotating high coverage separated contigs for $F3"

done

#------------------------------------------------------------------------------

for F4 in $(cat $list); do
mv results/09_high_coverage_contigs/tmp/$F4/*_tmp/*.gbk.txt.tmp results/09_high_coverage_contigs/tmp/$F4/
rm -rf results/09_high_coverage_contigs/tmp/$F4/*_tmp
done

#------------------------------------------------------------------------------

for F2 in $(cat $list); do
cd results/09_high_coverage_contigs/tmp/$F2

for F3 in $(cat list3.tmp); do
V1=$(grep $F3 *.csv.2xHCv.csv)
V2=$(cat *.HCv.contig.$F3.gbk.txt.tmp)
echo $V1 $V2 >> HCv.annotation.tmp
done

cd ..
cd ..
cd ..
cd ..

done

#------------------------------------------------------------------------------

for F1 in $(cat $list); do

awk -f /home/swapnil/pipeline/tools/vlookup.awk results/09_high_coverage_contigs/tmp/$F1/HCv.annotation.tmp results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.2xHCv.csv > results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.2xHCv.csv.labelled.tab
sed -i 's/==/ /g' results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.2xHCv.csv.labelled.tab

awk -f /home/swapnil/pipeline/tools/vlookup.awk results/09_high_coverage_contigs/tmp/$F1/HCv.annotation.tmp results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv > results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.labelled.tab
sed -i 's/==/ /g' results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.labelled.tab
cp results/09_high_coverage_contigs/tmp/$F1/$F1.10-500-filtered-contigs.csv.labelled.tab results/09_high_coverage_contigs/tmp/$F1.10-500-filtered-contigs.csv.labelled.tab

done

echo "finished.. separating HCv contigs and annotating"

## FINISH: separate high coverage contigs, annotate and write in csv
###############################################################################

echo "completed.. step-09 high coverage contigs analysis ----------------------"

###############################################################################


##discaraded code, dont know if useful

#------------------------------------------------------------------------------
## Writing the preliminary analysis results

#echo "running... writing the preliminary analysis results"

#for F1 in $(cat $list); do

#sed -i 's/, /-/g' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv

#------------------------------------------------------------------------------

## if the file doesnt contain translation column, add!

#if grep -q translation results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv; then
#:
#else
#awk '{print $0"SUFFIX"}' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.csv > results/09_high_coverage_contigs/tmp/$F1/tmp.csv
#sed 's/SUFFIX/\,"translation"/g' results/09_high_coverage_contigs/tmp/$F1/tmp.csv > results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.2.csv
#rm tmp.csv
#fi

#------------------------------------------------------------------------------

## separating each column

#awk -F',' 'NR==1{for (i=1;i<=NF;i++) a[i]=$i; next} {for (i=1;i<=NF;i++) {print $i > a[i]".'$F1'.tmp"}}' results/09_high_coverage_contigs/tmp/$F1/$F1.HCv.2.csv ; rename 's/"//' *.$F1.tmp ; rename 's/"././' *.$F1.tmp
#paste Start.$F1.tmp End.$F1.tmp Length.$F1.tmp locus_tag.$F1.tmp product.$F1.tmp translation.$F1.tmp >> results/09_high_coverage_contigs/all.HCv.csv

#rm *.$F1.tmp
#done
