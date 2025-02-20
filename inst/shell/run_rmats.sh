#!/usr/bash
# Parameters:
# $1:work directory
# $2:bam file directory
# $3:gtf file
# $4:paired or single
# $5:read length
# $6:core
# $7:path to rmats.py
# $8:path to python
set -e
workpath=$1
bampath=$2
gtf=$3
pair=$4
readlength=$5
core=$6
rMATS_path=$7
python_path=$8
mkdir -p $workpath

ls $bampath/*bam > $workpath/rMats.list
sed -i ':t;N;s/\n/,/;b t' $workpath/rMats.list

conda activate SCSES
python ${rMATS_path} --b1 $workpath/rMats.list --gtf $gtf -t $pair --variable-read-length \
        --readLength $readlength --nthread $core --od $workpath/output --tmp $workpath/tmp
rm -r $workpath/tmp
rm -r $workpath/output/tmp
rm $workpath/output/fromGTF* $workpath/output/JC*
