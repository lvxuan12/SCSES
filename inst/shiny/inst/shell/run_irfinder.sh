#!/usr/bash
# Parameters:
# $1:work directory
# $2:genome file
# $3:gtf file
# $4:bam file directory
# $5:core
# $6:readlength
# $7:path to IRFinder
# $8:path to samtools
# $9:path to star
set -e
workpath=$1
ref=$2
gtf=$3
bampath=$4
core=$5
readlength=$6
IRFinder_path=$7
samtools_path=$8
star_path=$9
#export LD_LIBRARY_PATH="/disk/lvxuan/lib:$LD_LIBRARY_PATH"

export PATH=$(dirname ${samtools_path}):$PATH

mkdir -p $workpath
cd $workpath

if [ -d "$workpath/REF" ]
then
    rm -r $workpath/REF
    echo "$workpath/REF should not yet exist! Remove $workpath/REF"
fi

$star_path --runMode genomeGenerate --runThreadN $core \
     --genomeDir $workpath/STAR_Reference \
     --genomeFastaFiles $ref \
     --sjdbGTFfile $gtf \
     --sjdbOverhang $((readlength-1))
star_ref=$workpath/STAR_Reference
$IRFinder_path BuildRefFromSTARRef -r $workpath/REF \
-x $star_ref \
-f $ref \
-g $gtf \
-t $core
rm -r $workpath/STAR_Reference

for bam in `ls $bampath/*bam`
do
    name=`echo $bam | awk -F'/' '{print $NF}' | sed 's/.bam//'`
    if [ -d "$workpath/$name" ]
    then
        rm -r $workpath/$name
        echo "$workpath/$name should not yet exist! Remove $workpath/$name"
    fi
    mkdir $workpath/$name
    ${samtools_path} sort -m 2G -n -@ $core -o $workpath/$name/${name}.bam $bam
    ${IRFinder_path} BAM -r $workpath/REF -d $workpath/$name $workpath/$name/${name}.bam && rm $workpath/$name/IRFinder-ChrCoverage.txt $workpath/$name/IRFinder-JuncCount.txt $workpath/$name/IRFinder-ROI.txt $workpath/$name/IRFinder-SpansPoint.txt $workpath/$name/${name}.bam
done
rm -r $workpath/REF


