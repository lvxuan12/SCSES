bam_file=$1
out_dir=$2
jar_path=$3
lib_path=$4
java_path=$5
core=$6
cell_info=$7
times=$8

bam_sc=$out_dir/bam/

n_row=$(wc -l $cell_info | awk '{print $1}')
n_row_split=$(echo $n_row| awk '{print int($1/"'$times'")}')
for i in $(seq 1 1 50)
do
    row_from=$(echo $i | awk '{print ($1-1)*"'${n_row_split}'"+1}')
    row_to=$(echo $i | awk '{print $1*"'${n_row_split}'"}')
    if [[ $i -eq 50 ]]
    then
        row_to=$n_row
    fi
    sed -n "${row_from},${row_to}p" $cell_info > ${cell_info}.${i}
    $java_path -Xmx100g -cp $lib_path/* -jar $jar_path $bam_file ${cell_info}.${i} $bam_sc $core
    echo "cell split ${row_from}-${row_to} done"
done



