#!/usr/bash
# Parameters:
# $1:work directory
# $2:genome file
# $3:gtf file
# $4:bam file directory
# $5:core
set -e
workpath=$1
ref=$2
gtf=$3
bampath=$4
core=$5
#export LD_LIBRARY_PATH="/disk/lvxuan/lib:$LD_LIBRARY_PATH"
mkdir -p $workpath
cd $workpath

if [ -d "$workpath/REF" ]
then
    rm -r $workpath/REF
    echo "$workpath/REF should not yet exist! Remove $workpath/REF"
fi
mkdir $workpath/REF
ln -s $ref $workpath/REF/genome.fa
ln -s $gtf $workpath/REF/transcripts.gtf 
IRFinder -m BuildRefProcess -r $workpath/REF
for bam in `ls $bampath/*bam`
do
    name=`echo $bam | awk -F'/' '{print $NF}' | sed 's/.bam//'`
    if [ -d "$workpath/$name" ]
    then
        rm -r $workpath/$name
        echo "$workpath/$name should not yet exist! Remove $workpath/REF"
    fi
    mkdir $workpath/$name
    samtools sort -m 2G -n -@ $core -o $workpath/$name/${name}.bam $bam
    IRFinder BAM -r $workpath/REF -d $workpath/$name $workpath/$name/${name}.bam && rm $workpath/$name/IRFinder-ChrCoverage.txt $workpath/$name/IRFinder-JuncCount.txt $workpath/$name/IRFinder-ROI.txt $workpath/$name/IRFinder-SpansPoint.txt $workpath/$name/${name}.bam
done
rm -r $workpath/REF


