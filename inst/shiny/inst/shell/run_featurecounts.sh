#!/usr/bash
# Parameters:
# $1:work directory
# $2:genome file
# $3:gtf file
# $4:bam file directory
# $5:core
# $6:output prefix
# $7:paired or single
# $8:path to featureCounts
workpath=$1
ref=$2
gtf=$3
bampath=$4
core=$5
dataSet=$6
pair=$7
featureCounts_path=$8
mkdir -p $workpath
bamfiles=$(ls $bampath/*bam)
if [[ $pair = "paired" ]]
then
    ${featureCounts_path} -T $core -p --countReadPairs -s 0 -M -O --fraction -Q 20 -t exon -g gene_name -a $gtf -o $workpath/${dataSet}_count.txt $bamfiles
else
    ${featureCounts_path} -T $core -s 0 -M -O --fraction -Q 20 -t exon -g gene_name -a $gtf -o $workpath/${dataSet}_count.txt $bamfiles
fi