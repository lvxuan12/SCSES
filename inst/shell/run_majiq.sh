#!/usr/bash
# Parameters:
# $1:work directory
# $2:bam file directory
# $3:gff file
# $4:core
# $5:MAJIQ ref name
# $6:Threshold on the minimum total number of reads for any junction
set -e
workpath=$1
bampath=$2
gff=$3
core=$4
genome_name=$5
junctionReads=$6
majiq_path=$7
voila_path=$8
license_file=$9
export MAJIQ_LICENSE_FILE=$license_file
mkdir -p $workpath
echo -e "[info]\nbamdirs=$bampath\ngenome=$genome_name\nstrandness=None\n\n[experiments]" > $workpath/majiq_build_config.ini
for f in `ls $bampath/*bam`
do
      name=`echo $f | awk -F'/' '{print $NF}' | sed 's/.bam//g'`
      echo "${name}=${name}" >> $workpath/majiq_build_config.ini
done
$majiq_path build $gff -c $workpath/majiq_build_config.ini -j$core -o $workpath/build --disable-ir --minreads $junctionReads
find $workpath/build -name *majiq > $workpath/majiq_path.txt
for i in `cat $workpath/majiq_path.txt`
do
      name=${i##*/}
      prefix=${name/\.majiq/}
      prefix2=$(echo $prefix | tr "-" "_")
      prefix3=$(echo $prefix2 | sed 's/0.//')
      echo $prefix3
      $majiq_path psi -o $workpath/psi -n $prefix3 $i
done
$voila_path modulize -d $workpath/modulize $workpath/build $workpath/psi -j$core --overwrite
rm -r $workpath/build $workpath/psi
