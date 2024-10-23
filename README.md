Menu
================

- [SCSES](#scses)
  - [Installation](#installation)
    - [Installation of dependencies and
      requirements](#installation-of-dependencies-and-requirements)
      - [1. python module](#1-python-module)
      - [2. Softwares](#2-softwares)
    - [Installation of SCSES](#installation-of-scses)
      - [Tips for some Installation
        errors](#tips-for-some-installation-errors)
  - [SCSES input](#scses-input)
  - [Getting started](#getting-started)
    - [Step0. Download and prepare test
      data](#step0-download-and-prepare-test-data)
    - [Step1. Read configure file](#step1-read-configure-file)
    - [Step2. Get gene expression](#step2-get-gene-expression)
      - [TPM matrix (for smart-seq2
        dataset)](#tpm-matrix-for-smart-seq2-dataset)
      - [Normalized UMI count matrix (for UMI
        dataset)](#normalized-umi-count-matrix-for-umi-dataset)
    - [Step3. Detect splicing events](#step3-detect-splicing-events)
      - [for smart-seq2 dataset](#for-smart-seq2-dataset)
      - [for UMI dataset](#for-umi-dataset)
    - [Step4. Quantify splicing events](#step4-quantify-splicing-events)
    - [Step5. Constructs similarity
      networks](#step5-constructs-similarity-networks)
      - [Parameters used in this step](#parameters-used-in-this-step)
    - [Step6. Imputation](#step6-imputation)
    - [Step7. Estimation](#step7-estimation)
      - [Fine-tune the model](#fine-tune-the-model)
    - [Step8. Cell Clustering](#step8-cell-clustering)
  - [Example](#example)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCSES

<!-- badges: start -->

<!-- badges: end -->

Single-cell Splicing Estimation based on Network Diffusion

## Installation

### Installation of dependencies and requirements

We recommend a new **conda** environment to install SCSES:

``` bash
conda create -n SCSES_test python=3.11
conda activate SCSES_test
```

To use SCSES, you will need to install R, Python, Matlab, and Java.

``` bash
conda install -c conda-forge r-base=4.3.1
```

#### 1. python module

``` bash
pip install pandas numpy scipy scikit-learn
pip install keras==2.15.0
pip install tensorflow==2.15.0.post1
```

#### 2. Softwares

To detect splicing events You will need to install rMATS, MAJIQ,
IRFinder. rMATS should be built in the same environment with SCSES (same
python). MAJIQ should be built in a new environment due to conflict of
python package version，

##### 2.1 [rMATS](https://github.com/Xinglab/rmats-turbo)

``` bash
wget https://github.com/Xinglab/rmats-turbo/releases/download/v4.3.0/rmats_turbo_v4_3_0.tar.gz
tar -zxvf rmats_turbo_v4_3_0.tar.gz
pip install Cython
./build_rmats
export PATH=/path/to/rmats_turbo_v4_3_0/:$PATH
```

##### 2.2 [MAJIQ](https://biociphers.bitbucket.io/majiq-docs/index.html)

``` bash
conda create -n MAJIQ python=3.11
conda activate MAJIQ
pip install git+https://bitbucket.org/biociphers/majiq_academic.git
export MAJIQ_LICENSE_FILE=/path/to/majiq_license_academic_official.lic
```

**NOTE:** MAJIQ will not function without providing the license file.

##### 2.3 [IRFinder](https://github.com/dgaolab/IRFinder)

``` bash
wget https://github.com/RitchieLabIGH/IRFinder/archive/refs/tags/v2.0.1.tar.gz
tar -zxvf v1.3.0.tar.gz
export PATH=/path/to/IRFinder-2.0.1/bin/:$PATH
```

**NOTE:** [STAR](https://github.com/alexdobin/STAR) is required to build
IRFinder reference! To run IRFinder correctly, you also need to install
STAR.

##### 2.4 [samtools](https://github.com/samtools/samtools)

### Installation of SCSES

Currently SCSES can only be installed from GitHub. To install SCSES,
type the following command in **R**:

``` r
options(download.file.method = "wget", times=100)
remotes::install_github("lvxuan12/SCSES")
```

#### Tips for some Installation errors

##### 1. cannot find fftw.h

``` bash
conda install conda-forge::fftw
```

##### 2. cannot find -lxml2

``` bash
conda install conda-forge::libxml2
```

##### 3. cannot find -lsz

``` bash
ln -s /usr/lib/x86_64-linux-gnu/libsz.so /path/to/miniconda/envs/SCSES_test/lib/libsz.so
```

## SCSES input

SCSES requires five inputs from the user

##### 1. bam files sorted by coordinate and index

##### 2. genome FASTA file and annotation GTF and GFF3 file

##### 3. configure file

A configure file is required to run SCSES. You can use
`createConfigshiny` command to generate a configure file:

``` r
library(SCSES)
createConfigshiny(host, port) 
```

After running this command, a interactive window will popup which allow
you to fill some parameters, such as Bam File Path, and Work Path. For a
detailed explanation of the configuration file, please refer to the
`config_anno.txt`.

Finally, you can click “Create Config” button and a json file will be
generated in the `work_path` you provided if successful.

##### 4. phast conservation file in bigWig format

For human and mouse, you could download it directly from UCSC browser:
[mm10.60way.phastCons.bw](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/),
[hg38.phastCons100way.bw](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/)
and
[hg19.100way.phastCons.bw](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/).

##### 5. RBP

Genes annotated as RBP are required to constructs similarity networks.
For human and mouse, you could read from `extdata/rbp`.

``` r
library(SCSES)
#> Loading required package: BiocManager
#> Bioconductor version '3.16' is out-of-date; the current release version '3.19'
#>   is available with R version '4.4'; see https://bioconductor.org/install
#> Warning: replacing previous import 'hdf5r::h5const' by 'rhdf5::h5const' when
#> loading 'SCSES'
#> Warning: replacing previous import 'hdf5r::h5version' by 'rhdf5::h5version'
#> when loading 'SCSES'
#> Warning: replacing previous import 'fs::path' by 'rtracklayer::path' when
#> loading 'SCSES'
#> Warning: replacing previous import 'hdf5r::values' by 'rtracklayer::values'
#> when loading 'SCSES'
# human
rbp = system.file("extdata/rbp/human_rbp.txt",package = "SCSES")
rbp <- readLines(rbp)
rbp[1:10]
#>  [1] "A1CF"       "AC004381.6" "ACIN1"      "ACO1"       "AKAP1"     
#>  [6] "ALKBH8"     "ALYREF"     "ANKHD1"     "ANKRD17"    "APTX"
```

## Getting started

### Step0. Download and prepare test data

The **test bam files** can be downloaded from
<https://doi.org/10.5281/zenodo.13943076>. The dataset includes three
cell lines (HCT1954, HepG2, and HL-60), with each cell type comprising
five cells. The cell identities can be found in `annotation.txt`.

The other input data required for SCSES for the test data can be
downloaded from <https://doi.org/10.5281/zenodo.13951695>, including
genome FASTA file, annotation GTF and GFF3 file, RBP file,STAR
reference, phast conservation file, and configure file.

Move test bam files to directory bam/. Move other input data to
directory refgenome/.

``` bash
ls /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam/
ls /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/
#> SRR11826368.bam
#> SRR11826368.bam.bai
#> SRR11826371.bam
#> SRR11826371.bam.bai
#> SRR11826409.bam
#> SRR11826409.bam.bai
#> SRR11826436.bam
#> SRR11826436.bam.bai
#> SRR11826437.bam
#> SRR11826437.bam.bai
#> SRR11826440.bam
#> SRR11826440.bam.bai
#> SRR11826455.bam
#> SRR11826455.bam.bai
#> SRR11826458.bam
#> SRR11826458.bam.bai
#> SRR11826464.bam
#> SRR11826464.bam.bai
#> SRR11826470.bam
#> SRR11826470.bam.bai
#> SRR1275148.bam
#> SRR1275148.bam.bai
#> SRR1275150.bam
#> SRR1275150.bam.bai
#> SRR1275227.bam
#> SRR1275227.bam.bai
#> SRR1275305.bam
#> SRR1275305.bam.bai
#> SRR1275317.bam
#> SRR1275317.bam.bai
#> STAR_Reference
#> annotation.txt
#> cell_line.json
#> hg19.100way.phastCons.bw
#> human_rbp.txt
#> test.fa
#> test.fa.fai
#> test.gff3
#> test.gtf
```

### Step1. Read configure file

``` r
library(SCSES)
#paras_file: path to configure file generated in the previous step
paras = readSCSESconfig(paras_file)
```

The `cell_line.json` file downloaded previously is an example
configuration file for test data. Users can modify this file or use
`createConfigshiny` function to create a new configuration file.

``` r
## Loading packages
library(SCSES)
#paras_file: path to configure file generated in the previous step
paras = readSCSESconfig('/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/cell_line.json')
names(paras)
#> [1] "DataSet" "Basic"   "Task"
print(paras$DataSet)
#> [1] "cell_line"
print(paras$Basic$bam_path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam"
print(paras$Basic$work_path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test"
```

### Step2. Get gene expression

#### TPM matrix (for smart-seq2 dataset)

The TPM matrix of gene expression can be obtained by different methods.
We use [featureCounts](https://github.com/ShiLab-Bioinformatics/subread)
to generate this matrix. You can use `getGeneExpression` to run
featureCounts, which will save featureCounts output to `work_path/expr/`
and `getEXPmatrix` to generate TPM matrix, which will save gene
expression count and TPM matrix to `work_path/rds/`.

``` r
featurecounts.path = getGeneExpression(paras) 
#> [1] "[2024-10-23 08:16:55] Detect gene expression: bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/shell/run_featurecounts.sh /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/expr/ /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/test.fa /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/test.gtf /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam 20 cell_line paired /disk/software/subread-2.0.6-source/bin/featureCounts >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/runfeatureCounts.log 2>&1"
#> [1] "[2024-10-23 08:17:14] Detect gene expression Finish."
rds.path = getEXPmatrix(paras)
print(rds.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds/"
list.files(rds.path)
#> [1] "count.rds" "event.rds" "psi.rds"   "rc.rds"    "TPM.rds"
```

``` r

tpm = readRDS(paste0(rds.path,'/TPM.rds'))
tpm[1:5,1:5]
#>            SRR11826368.bam SRR11826371.bam SRR11826409.bam SRR11826436.bam
#> DDX11L1           0.000000        0.000000        0.000000               0
#> WASH7P            2.215639        2.618154        4.034029               0
#> MIR1302-11        0.000000        0.000000        0.000000               0
#> FAM138A           0.000000        0.000000        0.000000               0
#> OR4G4P            0.000000        0.000000        0.000000               0
#>            SRR11826437.bam
#> DDX11L1                  0
#> WASH7P                   0
#> MIR1302-11               0
#> FAM138A                  0
#> OR4G4P                   0
```

#### Normalized UMI count matrix (for UMI dataset)

You can use `get10XEXPmatrix` to generate Normalized UMI count matrix
from 10X CellRanger hdf5 file, which will save normalized UMI count to
`work_path/rds/`.

``` r
rds.path = get10XEXPmatrix(paras,expr_path,sample_name)
```

### Step3. Detect splicing events

To define a global set of all splicing events, SCSES firstly merges all
bam files from every single cell to construct a pseudo-bulk bam file,
and identifies all types of splicing events by conventional algorithms.

#### for smart-seq2 dataset

``` r
pseudobulk.path = createPseudobulk(paras)
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data/"
#> [1] "[2024-10-23 08:17:15] Creating Pseudobulk directory..."
#> Warning in dir.create(path = pseudobulk.path, recursive = T):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data' already exists
#> [1] "Pseudobulk bam file exists."
#> [1] "[2024-10-23 08:17:15] Merge Bam Files: /disk/software/samtools/bin/samtools merge -f -@ 20 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data//all.bam /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam/*.bam --no-PG"
#> [1] "[2024-10-23 08:17:55] Merge Bam Files Finish."
#> [1] "[2024-10-23 08:17:55] Bam File Index: /disk/software/samtools/bin/samtools index -@ 20 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data//all.bam"
#> [1] "[2024-10-23 08:18:01] Bam File Index Finish."
print(pseudobulk.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data/"
list.files(pseudobulk.path)
#> [1] "all.bam"     "all.bam.bai"

#if you meet:
#irfinder: error while loading shared libraries: libboost_iostreams.so.1.71.0: cannot open shared object file: No such file or directory
old.ld=Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH = paste0("/disk/lvxuan/lib:", old.ld))
# you can Provide STAR reference to speed up the function
event.path = detectEvents(paras,star_ref_path="/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/STAR_Reference")
#> [1] "Checking cells..."
#> [1] "Checking events..."
#> [1] "event_type=SE;RI;A3SS;A5SS;MXE  checked"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/events/"
#> [1] "[2024-10-23 08:18:01] Creating events directory..."
#> Warning in dir.create(path = work_path):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/events' already
#> exists
#> [1] "[2024-10-23 08:18:01] Generating SE event id"
#> arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
#> [1] "[2024-10-23 08:18:02] Generating SE event id Finish."
#> [1] "[2024-10-23 08:18:02] Generating RI event id"
#> arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
#> [1] "[2024-10-23 08:18:16] Generating RI event id Finish."
#> [1] "[2024-10-23 08:18:16] Generating A3SS event id"
#> [1] "[2024-10-23 08:18:17] Generating A3SS event id Finish."
#> [1] "[2024-10-23 08:18:17] Generating A5SS event id"
#> [1] "[2024-10-23 08:18:17] Generating A5SS event id Finish."
#> [1] "[2024-10-23 08:18:17] Generating MXE event id"
#> [1] "[2024-10-23 08:18:17] Generating MXE event id Finish."
#> [1] "Total splicing event: SE=5764 RI=582 A3SS=237 A5SS=231 MXE=349"
print(event.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/events/"
list.files(event.path)
#>  [1] "A3SS.txt"            "A5SS.txt"            "IRFinder"           
#>  [4] "java_getRC_A3SS.log" "java_getRC_A5SS.log" "java_getRC_MXE.log" 
#>  [7] "java_getRC_RI.log"   "java_getRC_SE.log"   "majiq"              
#> [10] "MXE.txt"             "RI.txt"              "rMats"              
#> [13] "SE.txt"
```

Different types of splicing events will be saved to `work_path/events/`,
separately.

``` r
se.event=readLines(paste0(event.path,'SE.txt'))
print(nrow(se.event))
#> NULL
print(se.event[1])
#> [1] "isoform1=exon:chr5:180219870-180220097:-@junction:chr5:180220098-180229679:-@exon:chr5:180229680-180229803:-|isoform2=junction:chr5:180220098-180222655:-@exon:chr5:180222656-180222875:-@junction:chr5:180222876-180229679:-|MGAT1|SE"
```

#### for UMI dataset

SCSES requires single cell bam files being saved in a directory. For
UMI-based dataset, and using CellRanger for data process, the function
`split10XBAM` can be used to get single cell bam files.

``` r
# CellRanger_path: directory to CellRanger output
# out_path: directory to save single cell bam
# core: the number of threads
splitbam.path = split10XBAM(CellRanger_path,out_path,core)
# path to single cell bam files should be added to bam_path in configure file

# pseudobulk.path = createPseudobulk(paras)
# It is not necessary to execute `createPseudobulk`, and the `possorted_genome_bam.bam`,`possorted_genome_bam.bam.bai` from `CellRanger_path` can be moved to `work_path/data/all.bam`.
event.path = detectEvents(paras)
```

### Step4. Quantify splicing events

According to splicing events detected in the previous step, SCSES then
quantify raw reads associated with these splicing events in each cell,
and construct the raw read count matrix and calculate the raw PSI
matrix.

The definition of PSI of different AS events: ![The definition of PSI of
different AS events:](png/PSI.png)

``` r
rawrc.path = getRawRC(paras)
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "Checking cells..."
#> [1] "15 cells are considered."
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value/"
#> [1] "Splicing event types: A3SS A5SS MXE RI SE"
#> [1] "[2024-10-23 08:18:17] Counting reads of A3SS events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value//A3SS_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 100%"
#> [1] "[2024-10-23 08:18:20] Counting reads of A3SS events Finish."
#> [1] "[2024-10-23 08:18:20] Counting reads of A5SS events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value//A5SS_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 100%"
#> [1] "[2024-10-23 08:18:23] Counting reads of A5SS events Finish."
#> [1] "[2024-10-23 08:18:23] Counting reads of MXE events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value//MXE_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 100%"
#> [1] "[2024-10-23 08:18:28] Counting reads of MXE events Finish."
#> [1] "[2024-10-23 08:18:28] Counting reads of RI events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value//RI_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress:100%"
#> [1] "[2024-10-23 08:18:37] Counting reads of RI events Finish."
#> [1] "[2024-10-23 08:18:37] Counting reads of SE events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value//SE_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 100%"
#> [1] "[2024-10-23 08:19:06] Counting reads of SE events Finish."
rawpsi.path = getRawPSI(paras)
#> [1] "Checking raw reads..."
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of A3SS events..."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of A3SS events Finish."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of A5SS events..."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of A5SS events Finish."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of MXE events..."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of MXE events Finish."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of RI events..."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of RI events Finish."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of SE events..."
#> [1] "[2024-10-23 08:19:06] Calculating PSI value of SE events Finish."
rawrds.path = mergeSplicingValue(paras)
processed.data.path = preprocessEvent(paras)
#> [1] "[2024-10-23 08:19:07] Processing raw data..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds/"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds_processed/"
#> [1] "Before filtering"
#> [1] "expr: 17103*15"
#> [1] "psi: 7149*16"
#> [1] "rc: 16049*16"
#> [1] "Cell ids of raw data is mismatch! The insection of cells will be used to further analysis."
#> [1] "After filtering"
#> [1] "expr: 5669*15"
#> [1] "psi: 2550*15"
#> [1] "rc: 6087*15"
#> [1] "[2024-10-23 08:19:08] Successfully processed data."
print(rawrds.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds/"
print(processed.data.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds_processed/"
```

Raw read count matrix, PSI matrix, and event annotation will be saved to
`work_path/rds/`.Then, data after quality control process will be saved
to `work_path/rds_processed/`, which will be used for **subsequent
calculations**.

``` r
## read count for each type of splicing events in each single cell
list.files(rawrc.path,pattern = "*_rjm")
#> [1] "A3SS_rjm" "A5SS_rjm" "MXE_rjm"  "RI_rjm"   "SE_rjm"
## data for each type of splicing event
rawrc.files = list.files(rawrc.path,pattern = "*_rc.rds")
print(rawrc.files)
#> [1] "A3SS_rc.rds" "A5SS_rc.rds" "MXE_rc.rds"  "RI_rc.rds"   "SE_rc.rds"
rawpsi.files = list.files(rawpsi.path,pattern = "*_psi.rds")
print(rawpsi.files)
#> [1] "A3SS_psi.rds" "A5SS_psi.rds" "MXE_psi.rds"  "RI_psi.rds"   "SE_psi.rds"
## merged data
psi = readRDS(paste0(rawrds.path,'/psi.rds'))
print(psi[1:3,1:3])
#>                                                                                                                                                                                                                                       SRR11826368.bam
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS        0.8888889
#> isoform1=exon:chr2:233599865-233599948:+@junction:chr2:233599949-233600472:+@exon:chr2:233600473-233601358:+|isoform2=exon:chr2:233599865-233599948:+@junction:chr2:233599949-233600462:+@exon:chr2:233600463-233601358:+|GIGYF2|A3SS       0.0000000
#> isoform1=exon:chr2:233651860-233652039:+@junction:chr2:233652040-233655425:+@exon:chr2:233655426-233655625:+|isoform2=exon:chr2:233651860-233652039:+@junction:chr2:233652040-233655407:+@exon:chr2:233655408-233655625:+|GIGYF2|A3SS       1.0000000
#>                                                                                                                                                                                                                                       SRR11826371.bam
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS        1.0000000
#> isoform1=exon:chr2:233599865-233599948:+@junction:chr2:233599949-233600472:+@exon:chr2:233600473-233601358:+|isoform2=exon:chr2:233599865-233599948:+@junction:chr2:233599949-233600462:+@exon:chr2:233600463-233601358:+|GIGYF2|A3SS       0.0000000
#> isoform1=exon:chr2:233651860-233652039:+@junction:chr2:233652040-233655425:+@exon:chr2:233655426-233655625:+|isoform2=exon:chr2:233651860-233652039:+@junction:chr2:233652040-233655407:+@exon:chr2:233655408-233655625:+|GIGYF2|A3SS       0.6666667
#>                                                                                                                                                                                                                                       SRR11826409.bam
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS        0.8888889
#> isoform1=exon:chr2:233599865-233599948:+@junction:chr2:233599949-233600472:+@exon:chr2:233600473-233601358:+|isoform2=exon:chr2:233599865-233599948:+@junction:chr2:233599949-233600462:+@exon:chr2:233600463-233601358:+|GIGYF2|A3SS       0.0000000
#> isoform1=exon:chr2:233651860-233652039:+@junction:chr2:233652040-233655425:+@exon:chr2:233655426-233655625:+|isoform2=exon:chr2:233651860-233652039:+@junction:chr2:233652040-233655407:+@exon:chr2:233655408-233655625:+|GIGYF2|A3SS       0.9230769
print(dim(psi))
#> [1] 7149   16

rc = readRDS(paste0(rawrds.path,'/rc.rds'))
print(rc[1:3,1:3])
#>                                     SRR11826368.bam SRR11826371.bam
#> junction:chr4:146019572-146025540:+               2               0
#> junction:chr2:233599949-233600472:+               0               0
#> junction:chr2:233652040-233655425:+               0               2
#>                                     SRR11826409.bam
#> junction:chr4:146019572-146025540:+               1
#> junction:chr2:233599949-233600472:+               0
#> junction:chr2:233652040-233655425:+               1
print(dim(rc))
#> [1] 16049    16

event = readRDS(paste0(rawrds.path,'/event.rds'))
print(event[1:3,])
#>                                                                                                                                                                                                                                   event
#> 1  isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS
#> 2 isoform1=exon:chr2:233599865-233599948:+@junction:chr2:233599949-233600472:+@exon:chr2:233600473-233601358:+|isoform2=exon:chr2:233599865-233599948:+@junction:chr2:233599949-233600462:+@exon:chr2:233600463-233601358:+|GIGYF2|A3SS
#> 3 isoform1=exon:chr2:233651860-233652039:+@junction:chr2:233652040-233655425:+@exon:chr2:233655426-233655625:+|isoform2=exon:chr2:233651860-233652039:+@junction:chr2:233652040-233655407:+@exon:chr2:233655408-233655625:+|GIGYF2|A3SS
#>                            exclusion1 exclusion2
#> 1 junction:chr4:146019572-146025540:+       <NA>
#> 2 junction:chr2:233599949-233600472:+       <NA>
#> 3 junction:chr2:233652040-233655425:+       <NA>
#>                            retention1 retention2 type
#> 1 junction:chr4:146019572-146025537:+       <NA> A3SS
#> 2 junction:chr2:233599949-233600462:+       <NA> A3SS
#> 3 junction:chr2:233652040-233655407:+       <NA> A3SS

## processed data
psi_processed = readRDS(paste0(processed.data.path,'/psi.rds'))
rc_processed = readRDS(paste0(processed.data.path,'/rc.rds'))
print(dim(psi_processed))
#> [1] 2550   15
print(dim(rc_processed))
#> [1] 6087   15
```

### Step5. Constructs similarity networks

To overcome the high dropout rate and limited read coverage of scRNA-seq
techniques, SCSES constructs cell similarity and event similarity
networks by K-nearest neighbor algorithm (KNN) to learn information from
similar cells/events.

``` r
cellnet.path = getCellSimilarity(paras)
#> [1] "[2024-10-23 08:19:08] Calculate cell similarity..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds_processed/"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/cell_similarity/"
#> [1] "feature_num=1000  checked"
#> [1] "cell_similarity_data=PSI;RC;EXP_RBP  checked"
#> [1] "distance_method=euclidean  checked"
#> [1] "alpha_cell=0.8  checked"
#> [1] "kcell_max=8  checked"
#> [1] "kcell_min=3  checked"
#> [1] "decay_cell=0.05  checked"
#> [1] "Checking data: PSI"
#> [1] "2550 events are used to calculate cell similarity"
#> [1] "Checking data: RC"
#> [1] "6087 reads are used to calculate cell similarity"
#> [1] "Checking data: EXP_RBP"
#> [1] "384 rbps are used to calculate cell similarity"
#> [1] "[2024-10-23 08:19:08] Computing cell similarity based on PSI"
#> [1] "Calculate similarity among 15 cells."
#> Delta: [0.17081888 0.02261298]
#> [1] "[2024-10-23 08:19:29] Computing cell similarity based on RC"
#> [1] "Calculate similarity among 15 cells."
#> Delta: [0.13512509 0.02119545]
#> [1] "[2024-10-23 08:19:30] Computing cell similarity based on EXP_RBP"
#> [1] "Calculate similarity among 15 cells."
#> [1] "The number of features is greater than the number of rows in the input data."
#> [1] "Total 384 features will be used"
#> Delta: [0.17068148 0.03534559]
#> [1] "[2024-10-23 08:19:30] Calculate cell similarity Finish."
eventnet.path = getEventSimilarity(paras)
#> [1] "[2024-10-23 08:19:30] Calculate events KNN..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds_processed/"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/event_similarity/"
#> [1] "alpha_event=0.8  checked"
#> [1] "kevent=5  checked"
#> [1] "decay_event=0.05  checked"
#> [1] "Checking data..."
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "[2024-10-23 08:19:54] Calculate events feature..."
#> [1] "[2024-10-23 08:19:54] step1 Creating BSgenome for hg19======="
#> Warning: previous export ''hg19'' is being replaced
#> Warning: package 'hg19' was built under R version 4.3.1
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
#>     table, tapply, union, unique, unsplit, which.max, which.min
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> [1] "[2024-10-23 08:19:55] step2 Extracting features ======="
#> [1] "[2024-10-23 08:19:55] Extracting A3SS features..."
#> [1] "[2024-10-23 08:19:55] Loading events..."
#> [1] "[2024-10-23 08:20:06] Parsing events region..."
#> [1] "[2024-10-23 08:20:13] Extracting length features."
#> [1] "[2024-10-23 08:20:13] Extracting motif features."
#> [1] "[2024-10-23 08:20:13] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-23 08:20:14] Extracting kmer features."
#> [1] "[2024-10-23 08:20:14] Extracting A Ratio features."
#> [1] "[2024-10-23 08:20:15] Saving Result"
#> [1] "[2024-10-23 08:20:15] Extracting A3SS features Finished"
#> [1] "[2024-10-23 08:20:15] Extracting A5SS features..."
#> [1] "[2024-10-23 08:20:15] Loading events..."
#> [1] "[2024-10-23 08:20:26] Parsing events region..."
#> [1] "[2024-10-23 08:20:33] Extracting length features."
#> [1] "[2024-10-23 08:20:34] Extracting motif features."
#> [1] "[2024-10-23 08:20:34] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-23 08:20:34] Extracting kmer features."
#> [1] "[2024-10-23 08:20:34] Extracting A Ratio features."
#> [1] "[2024-10-23 08:20:36] Saving Result"
#> [1] "[2024-10-23 08:20:36] Extracting A5SS features Finished"
#> [1] "[2024-10-23 08:20:36] Extracting MXE features..."
#> [1] "[2024-10-23 08:20:36] Loading events..."
#> [1] "[2024-10-23 08:20:46] Parsing events region..."
#> [1] "[2024-10-23 08:20:55] Extracting length features."
#> [1] "[2024-10-23 08:20:55] Extracting motif features."
#> [1] "[2024-10-23 08:20:55] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-23 08:20:55] Extracting kmer features."
#> [1] "[2024-10-23 08:20:55] Extracting A Ratio features."
#> [1] "[2024-10-23 08:20:55] Saving Result"
#> [1] "[2024-10-23 08:20:55] Extracting MXE features Finished"
#> [1] "[2024-10-23 08:20:55] Extracting RI features..."
#> [1] "[2024-10-23 08:20:55] Loading events..."
#> [1] "[2024-10-23 08:21:06] Parsing events region..."
#> [1] "[2024-10-23 08:21:13] Extracting length features."
#> [1] "[2024-10-23 08:21:13] Extracting motif features."
#> [1] "[2024-10-23 08:21:13] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-23 08:21:14] Extracting kmer features."
#> [1] "[2024-10-23 08:21:14] Extracting A Ratio features."
#> [1] "[2024-10-23 08:21:15] Saving Result"
#> [1] "[2024-10-23 08:21:15] Extracting RI features Finished"
#> [1] "[2024-10-23 08:21:15] Extracting SE features..."
#> [1] "[2024-10-23 08:21:15] Loading events..."
#> [1] "[2024-10-23 08:21:26] Parsing events region..."
#> [1] "[2024-10-23 08:21:37] Extracting length features."
#> [1] "[2024-10-23 08:21:40] Extracting motif features."
#> [1] "[2024-10-23 08:21:40] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-23 08:21:49] Extracting kmer features."
#> [1] "[2024-10-23 08:21:49] Extracting A Ratio features."
#> [1] "[2024-10-23 08:21:49] Saving Result"
#> [1] "[2024-10-23 08:21:50] Extracting SE features Finished"
#> [1] "[2024-10-23 08:21:50] step3 Combining events feature ======="
#> [1] "[2024-10-23 08:21:50] Parsing A3SS features..."
#> [1] "[2024-10-23 08:21:50] Parsing A3SS features Finished"
#> [1] "[2024-10-23 08:21:50] Parsing A5SS features..."
#> [1] "[2024-10-23 08:21:50] Parsing A5SS features Finished"
#> [1] "[2024-10-23 08:21:50] Parsing MXE features..."
#> [1] "[2024-10-23 08:21:51] Parsing MXE features Finished"
#> [1] "[2024-10-23 08:21:51] Parsing RI features..."
#> [1] "[2024-10-23 08:21:51] Parsing RI features Finished"
#> [1] "[2024-10-23 08:21:51] Parsing SE features..."
#> [1] "[2024-10-23 08:21:53] Parsing SE features Finished"
#> [1] "[2024-10-23 08:21:53] step4 Encoding events feature ======="
#> [1] "[2024-10-23 08:21:53] A3SS event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02321/1 [==============================] - 2s 2s/step - loss: 1.0232 - val_loss: 1.0159
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01591/1 [==============================] - 0s 56ms/step - loss: 1.0159 - val_loss: 1.0131
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01311/1 [==============================] - 0s 69ms/step - loss: 1.0131 - val_loss: 1.0113
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01131/1 [==============================] - 0s 55ms/step - loss: 1.0113 - val_loss: 1.0095
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00951/1 [==============================] - 0s 63ms/step - loss: 1.0095 - val_loss: 1.0076
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00761/1 [==============================] - 0s 56ms/step - loss: 1.0076 - val_loss: 1.0051
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00511/1 [==============================] - 0s 52ms/step - loss: 1.0051 - val_loss: 1.0019
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00191/1 [==============================] - 0s 51ms/step - loss: 1.0019 - val_loss: 0.9975
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99751/1 [==============================] - 0s 53ms/step - loss: 0.9975 - val_loss: 0.9917
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99171/1 [==============================] - 0s 76ms/step - loss: 0.9917 - val_loss: 0.9842
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98421/1 [==============================] - 0s 60ms/step - loss: 0.9842 - val_loss: 0.9761
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97611/1 [==============================] - 0s 52ms/step - loss: 0.9761 - val_loss: 0.9686
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96861/1 [==============================] - 0s 72ms/step - loss: 0.9686 - val_loss: 0.9608
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96081/1 [==============================] - 0s 68ms/step - loss: 0.9608 - val_loss: 0.9514
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95141/1 [==============================] - 0s 57ms/step - loss: 0.9514 - val_loss: 0.9418
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94181/1 [==============================] - 0s 70ms/step - loss: 0.9418 - val_loss: 0.9321
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93211/1 [==============================] - 0s 63ms/step - loss: 0.9321 - val_loss: 0.9214
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92141/1 [==============================] - 0s 51ms/step - loss: 0.9214 - val_loss: 0.9109
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91091/1 [==============================] - 0s 74ms/step - loss: 0.9109 - val_loss: 0.9010
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90101/1 [==============================] - 0s 59ms/step - loss: 0.9010 - val_loss: 0.8908
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89081/1 [==============================] - 0s 72ms/step - loss: 0.8908 - val_loss: 0.8803
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88031/1 [==============================] - 0s 65ms/step - loss: 0.8803 - val_loss: 0.8694
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86941/1 [==============================] - 0s 54ms/step - loss: 0.8694 - val_loss: 0.8590
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85901/1 [==============================] - 0s 72ms/step - loss: 0.8590 - val_loss: 0.8487
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84871/1 [==============================] - 0s 70ms/step - loss: 0.8487 - val_loss: 0.8383
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83831/1 [==============================] - 0s 60ms/step - loss: 0.8383 - val_loss: 0.8282
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82821/1 [==============================] - 0s 72ms/step - loss: 0.8282 - val_loss: 0.8181
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81811/1 [==============================] - 0s 71ms/step - loss: 0.8181 - val_loss: 0.8082
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80821/1 [==============================] - 0s 74ms/step - loss: 0.8082 - val_loss: 0.7983
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79831/1 [==============================] - 0s 71ms/step - loss: 0.7983 - val_loss: 0.7884
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78841/1 [==============================] - 0s 56ms/step - loss: 0.7884 - val_loss: 0.7784
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77841/1 [==============================] - 0s 54ms/step - loss: 0.7784 - val_loss: 0.7684
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76841/1 [==============================] - 0s 73ms/step - loss: 0.7684 - val_loss: 0.7583
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75831/1 [==============================] - 0s 70ms/step - loss: 0.7583 - val_loss: 0.7484
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74841/1 [==============================] - 0s 56ms/step - loss: 0.7484 - val_loss: 0.7385
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73851/1 [==============================] - 0s 53ms/step - loss: 0.7385 - val_loss: 0.7285
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72851/1 [==============================] - 0s 53ms/step - loss: 0.7285 - val_loss: 0.7186
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71861/1 [==============================] - 0s 53ms/step - loss: 0.7186 - val_loss: 0.7088
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70881/1 [==============================] - 0s 69ms/step - loss: 0.7088 - val_loss: 0.6989
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69891/1 [==============================] - 0s 52ms/step - loss: 0.6989 - val_loss: 0.6892
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68921/1 [==============================] - 0s 61ms/step - loss: 0.6892 - val_loss: 0.6795
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67951/1 [==============================] - 0s 65ms/step - loss: 0.6795 - val_loss: 0.6699
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66991/1 [==============================] - 0s 67ms/step - loss: 0.6699 - val_loss: 0.6603
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66031/1 [==============================] - 0s 71ms/step - loss: 0.6603 - val_loss: 0.6508
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65081/1 [==============================] - 0s 65ms/step - loss: 0.6508 - val_loss: 0.6414
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64141/1 [==============================] - 0s 54ms/step - loss: 0.6414 - val_loss: 0.6320
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63201/1 [==============================] - 0s 53ms/step - loss: 0.6320 - val_loss: 0.6227
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62271/1 [==============================] - 0s 70ms/step - loss: 0.6227 - val_loss: 0.6135
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61351/1 [==============================] - 0s 57ms/step - loss: 0.6135 - val_loss: 0.6044
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60441/1 [==============================] - 0s 65ms/step - loss: 0.6044 - val_loss: 0.5954
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.59541/1 [==============================] - 0s 69ms/step - loss: 0.5954 - val_loss: 0.5864
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.58641/1 [==============================] - 0s 56ms/step - loss: 0.5864 - val_loss: 0.5775
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.57751/1 [==============================] - 0s 73ms/step - loss: 0.5775 - val_loss: 0.5687
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.56871/1 [==============================] - 0s 53ms/step - loss: 0.5687 - val_loss: 0.5599
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.55991/1 [==============================] - 0s 74ms/step - loss: 0.5599 - val_loss: 0.5512
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.55121/1 [==============================] - 0s 55ms/step - loss: 0.5512 - val_loss: 0.5426
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.54261/1 [==============================] - 0s 76ms/step - loss: 0.5426 - val_loss: 0.5339
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.53391/1 [==============================] - 0s 61ms/step - loss: 0.5339 - val_loss: 0.5254
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.52541/1 [==============================] - 0s 69ms/step - loss: 0.5254 - val_loss: 0.5169
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.51691/1 [==============================] - 0s 52ms/step - loss: 0.5169 - val_loss: 0.5085
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.50851/1 [==============================] - 0s 76ms/step - loss: 0.5085 - val_loss: 0.5001
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.50011/1 [==============================] - 0s 76ms/step - loss: 0.5001 - val_loss: 0.4918
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.49181/1 [==============================] - 0s 74ms/step - loss: 0.4918 - val_loss: 0.4835
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.48351/1 [==============================] - 0s 63ms/step - loss: 0.4835 - val_loss: 0.4753
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.47531/1 [==============================] - 0s 70ms/step - loss: 0.4753 - val_loss: 0.4672
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.46721/1 [==============================] - 0s 56ms/step - loss: 0.4672 - val_loss: 0.4591
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.45911/1 [==============================] - 0s 64ms/step - loss: 0.4591 - val_loss: 0.4511
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.45111/1 [==============================] - 0s 55ms/step - loss: 0.4511 - val_loss: 0.4432
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.44321/1 [==============================] - 0s 60ms/step - loss: 0.4432 - val_loss: 0.4353
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.43531/1 [==============================] - 0s 69ms/step - loss: 0.4353 - val_loss: 0.4274
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.42741/1 [==============================] - 0s 75ms/step - loss: 0.4274 - val_loss: 0.4197
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.41971/1 [==============================] - 0s 62ms/step - loss: 0.4197 - val_loss: 0.4121
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.41211/1 [==============================] - 0s 70ms/step - loss: 0.4121 - val_loss: 0.4045
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.40451/1 [==============================] - 0s 62ms/step - loss: 0.4045 - val_loss: 0.3970
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.39701/1 [==============================] - 0s 50ms/step - loss: 0.3970 - val_loss: 0.3896
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.38961/1 [==============================] - 0s 55ms/step - loss: 0.3896 - val_loss: 0.3823
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.38231/1 [==============================] - 0s 63ms/step - loss: 0.3823 - val_loss: 0.3750
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.37501/1 [==============================] - 0s 58ms/step - loss: 0.3750 - val_loss: 0.3680
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36801/1 [==============================] - 0s 59ms/step - loss: 0.3680 - val_loss: 0.3610
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36101/1 [==============================] - 0s 69ms/step - loss: 0.3610 - val_loss: 0.3540
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.35401/1 [==============================] - 0s 66ms/step - loss: 0.3540 - val_loss: 0.3473
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.34731/1 [==============================] - 0s 55ms/step - loss: 0.3473 - val_loss: 0.3406
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.34061/1 [==============================] - 0s 66ms/step - loss: 0.3406 - val_loss: 0.3340
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.33401/1 [==============================] - 0s 69ms/step - loss: 0.3340 - val_loss: 0.3276
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.32761/1 [==============================] - 0s 107ms/step - loss: 0.3276 - val_loss: 0.3213
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.32131/1 [==============================] - 0s 69ms/step - loss: 0.3213 - val_loss: 0.3151
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.31511/1 [==============================] - 0s 85ms/step - loss: 0.3151 - val_loss: 0.3090
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30901/1 [==============================] - 0s 69ms/step - loss: 0.3090 - val_loss: 0.3031
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30311/1 [==============================] - 0s 98ms/step - loss: 0.3031 - val_loss: 0.2973
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29731/1 [==============================] - 0s 62ms/step - loss: 0.2973 - val_loss: 0.2917
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29171/1 [==============================] - 0s 70ms/step - loss: 0.2917 - val_loss: 0.2861
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.28611/1 [==============================] - 0s 61ms/step - loss: 0.2861 - val_loss: 0.2807
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.28071/1 [==============================] - 0s 74ms/step - loss: 0.2807 - val_loss: 0.2753
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27531/1 [==============================] - 0s 70ms/step - loss: 0.2753 - val_loss: 0.2700
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27001/1 [==============================] - 0s 71ms/step - loss: 0.2700 - val_loss: 0.2649
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.26491/1 [==============================] - 0s 54ms/step - loss: 0.2649 - val_loss: 0.2599
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25991/1 [==============================] - 0s 71ms/step - loss: 0.2599 - val_loss: 0.2550
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25501/1 [==============================] - 0s 75ms/step - loss: 0.2550 - val_loss: 0.2502
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25021/1 [==============================] - 0s 60ms/step - loss: 0.2502 - val_loss: 0.2456
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.24561/1 [==============================] - 0s 70ms/step - loss: 0.2456 - val_loss: 0.2410
#> 1/4 [======>.......................] - ETA: 0s4/4 [==============================] - 0s 1ms/step
#> [1] "[2024-10-23 08:22:03] A3SS event encoding Finish."
#> [1] "[2024-10-23 08:22:03] A5SS event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02531/1 [==============================] - 1s 1s/step - loss: 1.0253 - val_loss: 1.0163
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01631/1 [==============================] - 0s 53ms/step - loss: 1.0163 - val_loss: 1.0134
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01341/1 [==============================] - 0s 50ms/step - loss: 1.0134 - val_loss: 1.0114
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01141/1 [==============================] - 0s 70ms/step - loss: 1.0114 - val_loss: 1.0096
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00961/1 [==============================] - 0s 52ms/step - loss: 1.0096 - val_loss: 1.0076
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00761/1 [==============================] - 0s 50ms/step - loss: 1.0076 - val_loss: 1.0052
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00521/1 [==============================] - 0s 51ms/step - loss: 1.0052 - val_loss: 1.0023
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00231/1 [==============================] - 0s 55ms/step - loss: 1.0023 - val_loss: 0.9985
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99851/1 [==============================] - 0s 69ms/step - loss: 0.9985 - val_loss: 0.9934
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99341/1 [==============================] - 0s 52ms/step - loss: 0.9934 - val_loss: 0.9869
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98691/1 [==============================] - 0s 68ms/step - loss: 0.9869 - val_loss: 0.9789
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97891/1 [==============================] - 0s 70ms/step - loss: 0.9789 - val_loss: 0.9697
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96971/1 [==============================] - 0s 68ms/step - loss: 0.9697 - val_loss: 0.9596
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95961/1 [==============================] - 0s 70ms/step - loss: 0.9596 - val_loss: 0.9478
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94781/1 [==============================] - 0s 56ms/step - loss: 0.9478 - val_loss: 0.9349
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93491/1 [==============================] - 0s 73ms/step - loss: 0.9349 - val_loss: 0.9221
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92211/1 [==============================] - 0s 70ms/step - loss: 0.9221 - val_loss: 0.9095
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90951/1 [==============================] - 0s 72ms/step - loss: 0.9095 - val_loss: 0.8970
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89701/1 [==============================] - 0s 67ms/step - loss: 0.8970 - val_loss: 0.8845
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88451/1 [==============================] - 0s 54ms/step - loss: 0.8845 - val_loss: 0.8718
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87181/1 [==============================] - 0s 50ms/step - loss: 0.8718 - val_loss: 0.8593
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85931/1 [==============================] - 0s 70ms/step - loss: 0.8593 - val_loss: 0.8476
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84761/1 [==============================] - 0s 52ms/step - loss: 0.8476 - val_loss: 0.8361
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83611/1 [==============================] - 0s 60ms/step - loss: 0.8361 - val_loss: 0.8246
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82461/1 [==============================] - 0s 56ms/step - loss: 0.8246 - val_loss: 0.8133
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81331/1 [==============================] - 0s 54ms/step - loss: 0.8133 - val_loss: 0.8022
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80221/1 [==============================] - 0s 80ms/step - loss: 0.8022 - val_loss: 0.7911
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79111/1 [==============================] - 0s 71ms/step - loss: 0.7911 - val_loss: 0.7800
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78001/1 [==============================] - 0s 53ms/step - loss: 0.7800 - val_loss: 0.7688
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76881/1 [==============================] - 0s 51ms/step - loss: 0.7688 - val_loss: 0.7577
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75771/1 [==============================] - 0s 53ms/step - loss: 0.7577 - val_loss: 0.7467
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74671/1 [==============================] - 0s 75ms/step - loss: 0.7467 - val_loss: 0.7357
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73571/1 [==============================] - 0s 55ms/step - loss: 0.7357 - val_loss: 0.7247
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72471/1 [==============================] - 0s 73ms/step - loss: 0.7247 - val_loss: 0.7138
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71381/1 [==============================] - 0s 68ms/step - loss: 0.7138 - val_loss: 0.7029
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70291/1 [==============================] - 0s 55ms/step - loss: 0.7029 - val_loss: 0.6920
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69201/1 [==============================] - 0s 51ms/step - loss: 0.6920 - val_loss: 0.6811
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68111/1 [==============================] - 0s 60ms/step - loss: 0.6811 - val_loss: 0.6704
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67041/1 [==============================] - 0s 74ms/step - loss: 0.6704 - val_loss: 0.6597
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65971/1 [==============================] - 0s 58ms/step - loss: 0.6597 - val_loss: 0.6491
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64911/1 [==============================] - 0s 50ms/step - loss: 0.6491 - val_loss: 0.6386
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63861/1 [==============================] - 0s 51ms/step - loss: 0.6386 - val_loss: 0.6282
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62821/1 [==============================] - 0s 50ms/step - loss: 0.6282 - val_loss: 0.6179
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61791/1 [==============================] - 0s 51ms/step - loss: 0.6179 - val_loss: 0.6077
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60771/1 [==============================] - 0s 54ms/step - loss: 0.6077 - val_loss: 0.5976
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.59761/1 [==============================] - 0s 50ms/step - loss: 0.5976 - val_loss: 0.5877
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.58771/1 [==============================] - 0s 52ms/step - loss: 0.5877 - val_loss: 0.5778
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.57781/1 [==============================] - 0s 51ms/step - loss: 0.5778 - val_loss: 0.5681
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.56811/1 [==============================] - 0s 52ms/step - loss: 0.5681 - val_loss: 0.5584
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.55841/1 [==============================] - 0s 57ms/step - loss: 0.5584 - val_loss: 0.5488
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.54881/1 [==============================] - 0s 76ms/step - loss: 0.5488 - val_loss: 0.5393
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.53931/1 [==============================] - 0s 56ms/step - loss: 0.5393 - val_loss: 0.5300
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.53001/1 [==============================] - 0s 72ms/step - loss: 0.5300 - val_loss: 0.5207
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.52071/1 [==============================] - 0s 55ms/step - loss: 0.5207 - val_loss: 0.5115
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.51151/1 [==============================] - 0s 56ms/step - loss: 0.5115 - val_loss: 0.5024
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.50241/1 [==============================] - 0s 80ms/step - loss: 0.5024 - val_loss: 0.4934
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.49341/1 [==============================] - 0s 72ms/step - loss: 0.4934 - val_loss: 0.4845
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.48451/1 [==============================] - 0s 59ms/step - loss: 0.4845 - val_loss: 0.4757
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.47571/1 [==============================] - 0s 53ms/step - loss: 0.4757 - val_loss: 0.4669
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.46691/1 [==============================] - 0s 53ms/step - loss: 0.4669 - val_loss: 0.4583
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.45831/1 [==============================] - 0s 76ms/step - loss: 0.4583 - val_loss: 0.4497
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.44971/1 [==============================] - 0s 74ms/step - loss: 0.4497 - val_loss: 0.4412
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.44121/1 [==============================] - 0s 71ms/step - loss: 0.4412 - val_loss: 0.4328
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.43281/1 [==============================] - 0s 76ms/step - loss: 0.4328 - val_loss: 0.4246
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.42461/1 [==============================] - 0s 73ms/step - loss: 0.4246 - val_loss: 0.4164
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.41641/1 [==============================] - 0s 72ms/step - loss: 0.4164 - val_loss: 0.4083
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.40831/1 [==============================] - 0s 56ms/step - loss: 0.4083 - val_loss: 0.4003
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.40031/1 [==============================] - 0s 73ms/step - loss: 0.4003 - val_loss: 0.3924
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.39241/1 [==============================] - 0s 73ms/step - loss: 0.3924 - val_loss: 0.3846
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.38461/1 [==============================] - 0s 75ms/step - loss: 0.3846 - val_loss: 0.3769
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.37691/1 [==============================] - 0s 69ms/step - loss: 0.3769 - val_loss: 0.3693
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36931/1 [==============================] - 0s 71ms/step - loss: 0.3693 - val_loss: 0.3619
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36191/1 [==============================] - 0s 69ms/step - loss: 0.3619 - val_loss: 0.3546
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.35461/1 [==============================] - 0s 69ms/step - loss: 0.3546 - val_loss: 0.3473
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.34731/1 [==============================] - 0s 53ms/step - loss: 0.3473 - val_loss: 0.3403
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.34031/1 [==============================] - 0s 71ms/step - loss: 0.3403 - val_loss: 0.3333
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.33331/1 [==============================] - 0s 58ms/step - loss: 0.3333 - val_loss: 0.3264
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.32641/1 [==============================] - 0s 51ms/step - loss: 0.3264 - val_loss: 0.3197
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.31971/1 [==============================] - 0s 50ms/step - loss: 0.3197 - val_loss: 0.3131
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.31311/1 [==============================] - 0s 51ms/step - loss: 0.3131 - val_loss: 0.3066
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30661/1 [==============================] - 0s 51ms/step - loss: 0.3066 - val_loss: 0.3002
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30021/1 [==============================] - 0s 71ms/step - loss: 0.3002 - val_loss: 0.2940
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29401/1 [==============================] - 0s 54ms/step - loss: 0.2940 - val_loss: 0.2879
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.28791/1 [==============================] - 0s 53ms/step - loss: 0.2879 - val_loss: 0.2820
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.28201/1 [==============================] - 0s 51ms/step - loss: 0.2820 - val_loss: 0.2761
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27611/1 [==============================] - 0s 52ms/step - loss: 0.2761 - val_loss: 0.2704
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27041/1 [==============================] - 0s 56ms/step - loss: 0.2704 - val_loss: 0.2648
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.26481/1 [==============================] - 0s 72ms/step - loss: 0.2648 - val_loss: 0.2593
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25931/1 [==============================] - 0s 56ms/step - loss: 0.2593 - val_loss: 0.2540
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25401/1 [==============================] - 0s 61ms/step - loss: 0.2540 - val_loss: 0.2488
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.24881/1 [==============================] - 0s 68ms/step - loss: 0.2488 - val_loss: 0.2436
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.24361/1 [==============================] - 0s 72ms/step - loss: 0.2436 - val_loss: 0.2386
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.23861/1 [==============================] - 0s 63ms/step - loss: 0.2386 - val_loss: 0.2338
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.23381/1 [==============================] - 0s 76ms/step - loss: 0.2338 - val_loss: 0.2290
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.22901/1 [==============================] - 0s 72ms/step - loss: 0.2290 - val_loss: 0.2243
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.22431/1 [==============================] - 0s 60ms/step - loss: 0.2243 - val_loss: 0.2197
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.21971/1 [==============================] - 0s 51ms/step - loss: 0.2197 - val_loss: 0.2153
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.21531/1 [==============================] - 0s 52ms/step - loss: 0.2153 - val_loss: 0.2109
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.21091/1 [==============================] - 0s 50ms/step - loss: 0.2109 - val_loss: 0.2067
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.20671/1 [==============================] - 0s 50ms/step - loss: 0.2067 - val_loss: 0.2025
#> 1/4 [======>.......................] - ETA: 0s4/4 [==============================] - 0s 1ms/step
#> [1] "[2024-10-23 08:22:11] A5SS event encoding Finish."
#> [1] "[2024-10-23 08:22:11] MXE event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99001/1 [==============================] - 1s 1s/step - loss: 0.9900 - val_loss: 0.9778
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97781/1 [==============================] - 0s 53ms/step - loss: 0.9778 - val_loss: 0.9738
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97381/1 [==============================] - 0s 68ms/step - loss: 0.9738 - val_loss: 0.9712
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97121/1 [==============================] - 0s 51ms/step - loss: 0.9712 - val_loss: 0.9683
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96831/1 [==============================] - 0s 69ms/step - loss: 0.9683 - val_loss: 0.9652
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96521/1 [==============================] - 0s 68ms/step - loss: 0.9652 - val_loss: 0.9614
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96141/1 [==============================] - 0s 53ms/step - loss: 0.9614 - val_loss: 0.9565
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95651/1 [==============================] - 0s 70ms/step - loss: 0.9565 - val_loss: 0.9502
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95021/1 [==============================] - 0s 62ms/step - loss: 0.9502 - val_loss: 0.9419
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94191/1 [==============================] - 0s 63ms/step - loss: 0.9419 - val_loss: 0.9311
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93111/1 [==============================] - 0s 54ms/step - loss: 0.9311 - val_loss: 0.9178
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91781/1 [==============================] - 0s 62ms/step - loss: 0.9178 - val_loss: 0.9026
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90261/1 [==============================] - 0s 72ms/step - loss: 0.9026 - val_loss: 0.8866
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88661/1 [==============================] - 0s 67ms/step - loss: 0.8866 - val_loss: 0.8705
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87051/1 [==============================] - 0s 54ms/step - loss: 0.8705 - val_loss: 0.8537
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85371/1 [==============================] - 0s 48ms/step - loss: 0.8537 - val_loss: 0.8348
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83481/1 [==============================] - 0s 50ms/step - loss: 0.8348 - val_loss: 0.8145
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81451/1 [==============================] - 0s 68ms/step - loss: 0.8145 - val_loss: 0.7938
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79381/1 [==============================] - 0s 53ms/step - loss: 0.7938 - val_loss: 0.7733
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77331/1 [==============================] - 0s 51ms/step - loss: 0.7733 - val_loss: 0.7524
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75241/1 [==============================] - 0s 54ms/step - loss: 0.7524 - val_loss: 0.7311
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73111/1 [==============================] - 0s 76ms/step - loss: 0.7311 - val_loss: 0.7097
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70971/1 [==============================] - 0s 69ms/step - loss: 0.7097 - val_loss: 0.6883
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68831/1 [==============================] - 0s 74ms/step - loss: 0.6883 - val_loss: 0.6671
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66711/1 [==============================] - 0s 60ms/step - loss: 0.6671 - val_loss: 0.6460
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64601/1 [==============================] - 0s 61ms/step - loss: 0.6460 - val_loss: 0.6251
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62511/1 [==============================] - 0s 65ms/step - loss: 0.6251 - val_loss: 0.6045
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60451/1 [==============================] - 0s 63ms/step - loss: 0.6045 - val_loss: 0.5837
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.58371/1 [==============================] - 0s 49ms/step - loss: 0.5837 - val_loss: 0.5631
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.56311/1 [==============================] - 0s 58ms/step - loss: 0.5631 - val_loss: 0.5429
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.54291/1 [==============================] - 0s 68ms/step - loss: 0.5429 - val_loss: 0.5232
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.52321/1 [==============================] - 0s 65ms/step - loss: 0.5232 - val_loss: 0.5038
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.50381/1 [==============================] - 0s 71ms/step - loss: 0.5038 - val_loss: 0.4847
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.48471/1 [==============================] - 0s 67ms/step - loss: 0.4847 - val_loss: 0.4662
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.46621/1 [==============================] - 0s 52ms/step - loss: 0.4662 - val_loss: 0.4482
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.44821/1 [==============================] - 0s 67ms/step - loss: 0.4482 - val_loss: 0.4307
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.43071/1 [==============================] - 0s 50ms/step - loss: 0.4307 - val_loss: 0.4138
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.41381/1 [==============================] - 0s 47ms/step - loss: 0.4138 - val_loss: 0.3973
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.39731/1 [==============================] - 0s 47ms/step - loss: 0.3973 - val_loss: 0.3816
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.38161/1 [==============================] - 0s 49ms/step - loss: 0.3816 - val_loss: 0.3664
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36641/1 [==============================] - 0s 64ms/step - loss: 0.3664 - val_loss: 0.3518
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.35181/1 [==============================] - 0s 49ms/step - loss: 0.3518 - val_loss: 0.3377
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.33771/1 [==============================] - 0s 47ms/step - loss: 0.3377 - val_loss: 0.3243
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.32431/1 [==============================] - 0s 49ms/step - loss: 0.3243 - val_loss: 0.3113
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.31131/1 [==============================] - 0s 66ms/step - loss: 0.3113 - val_loss: 0.2988
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29881/1 [==============================] - 0s 74ms/step - loss: 0.2988 - val_loss: 0.2869
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.28691/1 [==============================] - 0s 67ms/step - loss: 0.2869 - val_loss: 0.2754
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27541/1 [==============================] - 0s 63ms/step - loss: 0.2754 - val_loss: 0.2644
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.26441/1 [==============================] - 0s 66ms/step - loss: 0.2644 - val_loss: 0.2537
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25371/1 [==============================] - 0s 53ms/step - loss: 0.2537 - val_loss: 0.2436
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.24361/1 [==============================] - 0s 62ms/step - loss: 0.2436 - val_loss: 0.2339
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.23391/1 [==============================] - 0s 51ms/step - loss: 0.2339 - val_loss: 0.2245
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.22451/1 [==============================] - 0s 65ms/step - loss: 0.2245 - val_loss: 0.2156
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.21561/1 [==============================] - 0s 59ms/step - loss: 0.2156 - val_loss: 0.2070
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.20701/1 [==============================] - 0s 66ms/step - loss: 0.2070 - val_loss: 0.1989
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.19891/1 [==============================] - 0s 53ms/step - loss: 0.1989 - val_loss: 0.1911
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.19111/1 [==============================] - 0s 66ms/step - loss: 0.1911 - val_loss: 0.1837
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.18371/1 [==============================] - 0s 59ms/step - loss: 0.1837 - val_loss: 0.1765
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.17651/1 [==============================] - 0s 49ms/step - loss: 0.1765 - val_loss: 0.1698
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.16981/1 [==============================] - 0s 59ms/step - loss: 0.1698 - val_loss: 0.1632
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.16321/1 [==============================] - 0s 62ms/step - loss: 0.1632 - val_loss: 0.1570
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.15701/1 [==============================] - 0s 53ms/step - loss: 0.1570 - val_loss: 0.1510
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.15101/1 [==============================] - 0s 47ms/step - loss: 0.1510 - val_loss: 0.1453
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.14531/1 [==============================] - 0s 47ms/step - loss: 0.1453 - val_loss: 0.1398
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.13981/1 [==============================] - 0s 51ms/step - loss: 0.1398 - val_loss: 0.1345
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.13451/1 [==============================] - 0s 46ms/step - loss: 0.1345 - val_loss: 0.1293
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.12931/1 [==============================] - 0s 60ms/step - loss: 0.1293 - val_loss: 0.1244
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.12441/1 [==============================] - 0s 56ms/step - loss: 0.1244 - val_loss: 0.1196
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.11961/1 [==============================] - 0s 54ms/step - loss: 0.1196 - val_loss: 0.1151
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.11511/1 [==============================] - 0s 68ms/step - loss: 0.1151 - val_loss: 0.1107
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.11071/1 [==============================] - 0s 68ms/step - loss: 0.1107 - val_loss: 0.1064
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.10641/1 [==============================] - 0s 66ms/step - loss: 0.1064 - val_loss: 0.1023
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.10231/1 [==============================] - 0s 53ms/step - loss: 0.1023 - val_loss: 0.0984
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.09841/1 [==============================] - 0s 68ms/step - loss: 0.0984 - val_loss: 0.0946
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.09461/1 [==============================] - 0s 53ms/step - loss: 0.0946 - val_loss: 0.0910
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.09101/1 [==============================] - 0s 60ms/step - loss: 0.0910 - val_loss: 0.0876
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.08761/1 [==============================] - 0s 55ms/step - loss: 0.0876 - val_loss: 0.0844
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.08441/1 [==============================] - 0s 69ms/step - loss: 0.0844 - val_loss: 0.0813
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.08131/1 [==============================] - 0s 52ms/step - loss: 0.0813 - val_loss: 0.0783
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.07831/1 [==============================] - 0s 47ms/step - loss: 0.0783 - val_loss: 0.0756
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.07561/1 [==============================] - 0s 54ms/step - loss: 0.0756 - val_loss: 0.0730
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.07301/1 [==============================] - 0s 65ms/step - loss: 0.0730 - val_loss: 0.0706
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.07061/1 [==============================] - 0s 50ms/step - loss: 0.0706 - val_loss: 0.0684
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06841/1 [==============================] - 0s 49ms/step - loss: 0.0684 - val_loss: 0.0663
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06631/1 [==============================] - 0s 68ms/step - loss: 0.0663 - val_loss: 0.0643
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06431/1 [==============================] - 0s 66ms/step - loss: 0.0643 - val_loss: 0.0625
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06251/1 [==============================] - 0s 64ms/step - loss: 0.0625 - val_loss: 0.0609
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06091/1 [==============================] - 0s 66ms/step - loss: 0.0609 - val_loss: 0.0594
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05941/1 [==============================] - 0s 52ms/step - loss: 0.0594 - val_loss: 0.0580
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05801/1 [==============================] - 0s 49ms/step - loss: 0.0580 - val_loss: 0.0567
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05671/1 [==============================] - 0s 73ms/step - loss: 0.0567 - val_loss: 0.0555
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05551/1 [==============================] - 0s 58ms/step - loss: 0.0555 - val_loss: 0.0544
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05441/1 [==============================] - 0s 47ms/step - loss: 0.0544 - val_loss: 0.0534
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05341/1 [==============================] - 0s 46ms/step - loss: 0.0534 - val_loss: 0.0525
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05251/1 [==============================] - 0s 47ms/step - loss: 0.0525 - val_loss: 0.0517
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05171/1 [==============================] - 0s 46ms/step - loss: 0.0517 - val_loss: 0.0509
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05091/1 [==============================] - 0s 48ms/step - loss: 0.0509 - val_loss: 0.0502
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05021/1 [==============================] - 0s 49ms/step - loss: 0.0502 - val_loss: 0.0495
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.04951/1 [==============================] - 0s 47ms/step - loss: 0.0495 - val_loss: 0.0488
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.04881/1 [==============================] - 0s 47ms/step - loss: 0.0488 - val_loss: 0.0482
#> 1/2 [==============>...............] - ETA: 0s2/2 [==============================] - 0s 1ms/step
#> [1] "[2024-10-23 08:22:18] MXE event encoding Finish."
#> [1] "[2024-10-23 08:22:18] RI event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02971/1 [==============================] - 1s 1s/step - loss: 1.0297 - val_loss: 1.0206
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02061/1 [==============================] - 0s 58ms/step - loss: 1.0206 - val_loss: 1.0169
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01691/1 [==============================] - 0s 56ms/step - loss: 1.0169 - val_loss: 1.0143
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01431/1 [==============================] - 0s 57ms/step - loss: 1.0143 - val_loss: 1.0109
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01091/1 [==============================] - 0s 56ms/step - loss: 1.0109 - val_loss: 1.0060
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00601/1 [==============================] - 0s 58ms/step - loss: 1.0060 - val_loss: 0.9990
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99901/1 [==============================] - 0s 56ms/step - loss: 0.9990 - val_loss: 0.9905
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99051/1 [==============================] - 0s 57ms/step - loss: 0.9905 - val_loss: 0.9821
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98211/1 [==============================] - 0s 64ms/step - loss: 0.9821 - val_loss: 0.9758
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97581/1 [==============================] - 0s 68ms/step - loss: 0.9758 - val_loss: 0.9680
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96801/1 [==============================] - 0s 72ms/step - loss: 0.9680 - val_loss: 0.9590
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95901/1 [==============================] - 0s 78ms/step - loss: 0.9590 - val_loss: 0.9520
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95201/1 [==============================] - 0s 67ms/step - loss: 0.9520 - val_loss: 0.9466
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94661/1 [==============================] - 0s 75ms/step - loss: 0.9466 - val_loss: 0.9413
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94131/1 [==============================] - 0s 68ms/step - loss: 0.9413 - val_loss: 0.9364
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93641/1 [==============================] - 0s 74ms/step - loss: 0.9364 - val_loss: 0.9318
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93181/1 [==============================] - 0s 65ms/step - loss: 0.9318 - val_loss: 0.9256
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92561/1 [==============================] - 0s 86ms/step - loss: 0.9256 - val_loss: 0.9186
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91861/1 [==============================] - 0s 73ms/step - loss: 0.9186 - val_loss: 0.9116
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91161/1 [==============================] - 0s 66ms/step - loss: 0.9116 - val_loss: 0.9047
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90471/1 [==============================] - 0s 84ms/step - loss: 0.9047 - val_loss: 0.9013
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90131/1 [==============================] - 0s 62ms/step - loss: 0.9013 - val_loss: 0.8980
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89801/1 [==============================] - 0s 75ms/step - loss: 0.8980 - val_loss: 0.8913
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89131/1 [==============================] - 0s 78ms/step - loss: 0.8913 - val_loss: 0.8847
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88471/1 [==============================] - 0s 74ms/step - loss: 0.8847 - val_loss: 0.8797
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87971/1 [==============================] - 0s 77ms/step - loss: 0.8797 - val_loss: 0.8760
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87601/1 [==============================] - 0s 62ms/step - loss: 0.8760 - val_loss: 0.8722
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87221/1 [==============================] - 0s 69ms/step - loss: 0.8722 - val_loss: 0.8674
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86741/1 [==============================] - 0s 68ms/step - loss: 0.8674 - val_loss: 0.8629
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86291/1 [==============================] - 0s 75ms/step - loss: 0.8629 - val_loss: 0.8593
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85931/1 [==============================] - 0s 57ms/step - loss: 0.8593 - val_loss: 0.8564
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85641/1 [==============================] - 0s 60ms/step - loss: 0.8564 - val_loss: 0.8527
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85271/1 [==============================] - 0s 75ms/step - loss: 0.8527 - val_loss: 0.8484
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84841/1 [==============================] - 0s 66ms/step - loss: 0.8484 - val_loss: 0.8449
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84491/1 [==============================] - 0s 68ms/step - loss: 0.8449 - val_loss: 0.8416
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84161/1 [==============================] - 0s 59ms/step - loss: 0.8416 - val_loss: 0.8380
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83801/1 [==============================] - 0s 79ms/step - loss: 0.8380 - val_loss: 0.8340
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83401/1 [==============================] - 0s 67ms/step - loss: 0.8340 - val_loss: 0.8303
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83031/1 [==============================] - 0s 70ms/step - loss: 0.8303 - val_loss: 0.8267
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82671/1 [==============================] - 0s 71ms/step - loss: 0.8267 - val_loss: 0.8228
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82281/1 [==============================] - 0s 54ms/step - loss: 0.8228 - val_loss: 0.8186
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81861/1 [==============================] - 0s 57ms/step - loss: 0.8186 - val_loss: 0.8146
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81461/1 [==============================] - 0s 67ms/step - loss: 0.8146 - val_loss: 0.8106
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81061/1 [==============================] - 0s 79ms/step - loss: 0.8106 - val_loss: 0.8064
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80641/1 [==============================] - 0s 76ms/step - loss: 0.8064 - val_loss: 0.8022
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80221/1 [==============================] - 0s 80ms/step - loss: 0.8022 - val_loss: 0.7980
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79801/1 [==============================] - 0s 75ms/step - loss: 0.7980 - val_loss: 0.7938
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79381/1 [==============================] - 0s 72ms/step - loss: 0.7938 - val_loss: 0.7895
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78951/1 [==============================] - 0s 64ms/step - loss: 0.7895 - val_loss: 0.7852
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78521/1 [==============================] - 0s 72ms/step - loss: 0.7852 - val_loss: 0.7810
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78101/1 [==============================] - 0s 74ms/step - loss: 0.7810 - val_loss: 0.7767
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77671/1 [==============================] - 0s 60ms/step - loss: 0.7767 - val_loss: 0.7726
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77261/1 [==============================] - 0s 75ms/step - loss: 0.7726 - val_loss: 0.7685
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76851/1 [==============================] - 0s 59ms/step - loss: 0.7685 - val_loss: 0.7645
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76451/1 [==============================] - 0s 59ms/step - loss: 0.7645 - val_loss: 0.7605
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76051/1 [==============================] - 0s 73ms/step - loss: 0.7605 - val_loss: 0.7566
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75661/1 [==============================] - 0s 57ms/step - loss: 0.7566 - val_loss: 0.7527
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75271/1 [==============================] - 0s 55ms/step - loss: 0.7527 - val_loss: 0.7488
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74881/1 [==============================] - 0s 71ms/step - loss: 0.7488 - val_loss: 0.7449
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74491/1 [==============================] - 0s 61ms/step - loss: 0.7449 - val_loss: 0.7411
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74111/1 [==============================] - 0s 55ms/step - loss: 0.7411 - val_loss: 0.7373
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73731/1 [==============================] - 0s 66ms/step - loss: 0.7373 - val_loss: 0.7336
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73361/1 [==============================] - 0s 57ms/step - loss: 0.7336 - val_loss: 0.7298
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72981/1 [==============================] - 0s 75ms/step - loss: 0.7298 - val_loss: 0.7261
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72611/1 [==============================] - 0s 73ms/step - loss: 0.7261 - val_loss: 0.7225
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72251/1 [==============================] - 0s 72ms/step - loss: 0.7225 - val_loss: 0.7188
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71881/1 [==============================] - 0s 71ms/step - loss: 0.7188 - val_loss: 0.7151
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71511/1 [==============================] - 0s 74ms/step - loss: 0.7151 - val_loss: 0.7115
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71151/1 [==============================] - 0s 75ms/step - loss: 0.7115 - val_loss: 0.7079
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70791/1 [==============================] - 0s 71ms/step - loss: 0.7079 - val_loss: 0.7043
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70431/1 [==============================] - 0s 59ms/step - loss: 0.7043 - val_loss: 0.7007
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70071/1 [==============================] - 0s 64ms/step - loss: 0.7007 - val_loss: 0.6972
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69721/1 [==============================] - 0s 74ms/step - loss: 0.6972 - val_loss: 0.6938
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69381/1 [==============================] - 0s 79ms/step - loss: 0.6938 - val_loss: 0.6904
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69041/1 [==============================] - 0s 72ms/step - loss: 0.6904 - val_loss: 0.6869
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68691/1 [==============================] - 0s 68ms/step - loss: 0.6869 - val_loss: 0.6832
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68321/1 [==============================] - 0s 73ms/step - loss: 0.6832 - val_loss: 0.6796
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67961/1 [==============================] - 0s 68ms/step - loss: 0.6796 - val_loss: 0.6762
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67621/1 [==============================] - 0s 62ms/step - loss: 0.6762 - val_loss: 0.6729
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67291/1 [==============================] - 0s 56ms/step - loss: 0.6729 - val_loss: 0.6696
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66961/1 [==============================] - 0s 54ms/step - loss: 0.6696 - val_loss: 0.6662
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66621/1 [==============================] - 0s 56ms/step - loss: 0.6662 - val_loss: 0.6627
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66271/1 [==============================] - 0s 74ms/step - loss: 0.6627 - val_loss: 0.6593
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65931/1 [==============================] - 0s 77ms/step - loss: 0.6593 - val_loss: 0.6559
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65591/1 [==============================] - 0s 73ms/step - loss: 0.6559 - val_loss: 0.6527
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65271/1 [==============================] - 0s 57ms/step - loss: 0.6527 - val_loss: 0.6494
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64941/1 [==============================] - 0s 55ms/step - loss: 0.6494 - val_loss: 0.6461
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64611/1 [==============================] - 0s 67ms/step - loss: 0.6461 - val_loss: 0.6429
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64291/1 [==============================] - 0s 67ms/step - loss: 0.6429 - val_loss: 0.6396
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63961/1 [==============================] - 0s 54ms/step - loss: 0.6396 - val_loss: 0.6363
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63631/1 [==============================] - 0s 74ms/step - loss: 0.6363 - val_loss: 0.6330
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63301/1 [==============================] - 0s 63ms/step - loss: 0.6330 - val_loss: 0.6296
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62961/1 [==============================] - 0s 77ms/step - loss: 0.6296 - val_loss: 0.6264
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62641/1 [==============================] - 0s 70ms/step - loss: 0.6264 - val_loss: 0.6233
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62331/1 [==============================] - 0s 72ms/step - loss: 0.6233 - val_loss: 0.6201
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62011/1 [==============================] - 0s 57ms/step - loss: 0.6201 - val_loss: 0.6170
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61701/1 [==============================] - 0s 73ms/step - loss: 0.6170 - val_loss: 0.6137
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61371/1 [==============================] - 0s 71ms/step - loss: 0.6137 - val_loss: 0.6104
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61041/1 [==============================] - 0s 56ms/step - loss: 0.6104 - val_loss: 0.6072
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60721/1 [==============================] - 0s 75ms/step - loss: 0.6072 - val_loss: 0.6041
#>  1/12 [=>............................] - ETA: 0s12/12 [==============================] - 0s 1ms/step
#> [1] "[2024-10-23 08:22:27] RI event encoding Finish."
#> [1] "[2024-10-23 08:22:27] SE event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.03361/1 [==============================] - 1s 1s/step - loss: 1.0336 - val_loss: 1.0230
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02301/1 [==============================] - 0s 72ms/step - loss: 1.0230 - val_loss: 1.0188
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01881/1 [==============================] - 0s 73ms/step - loss: 1.0188 - val_loss: 1.0166
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01661/1 [==============================] - 0s 77ms/step - loss: 1.0166 - val_loss: 1.0152
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01521/1 [==============================] - 0s 81ms/step - loss: 1.0152 - val_loss: 1.0141
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01411/1 [==============================] - 0s 70ms/step - loss: 1.0141 - val_loss: 1.0130
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01301/1 [==============================] - 0s 76ms/step - loss: 1.0130 - val_loss: 1.0119
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01191/1 [==============================] - 0s 67ms/step - loss: 1.0119 - val_loss: 1.0106
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01061/1 [==============================] - 0s 77ms/step - loss: 1.0106 - val_loss: 1.0089
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00891/1 [==============================] - 0s 96ms/step - loss: 1.0089 - val_loss: 1.0067
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00671/1 [==============================] - 0s 76ms/step - loss: 1.0067 - val_loss: 1.0039
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00391/1 [==============================] - 0s 68ms/step - loss: 1.0039 - val_loss: 1.0005
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00051/1 [==============================] - 0s 69ms/step - loss: 1.0005 - val_loss: 0.9969
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99691/1 [==============================] - 0s 91ms/step - loss: 0.9969 - val_loss: 0.9940
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99401/1 [==============================] - 0s 73ms/step - loss: 0.9940 - val_loss: 0.9914
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99141/1 [==============================] - 0s 81ms/step - loss: 0.9914 - val_loss: 0.9879
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98791/1 [==============================] - 0s 77ms/step - loss: 0.9879 - val_loss: 0.9836
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98361/1 [==============================] - 0s 71ms/step - loss: 0.9836 - val_loss: 0.9795
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97951/1 [==============================] - 0s 70ms/step - loss: 0.9795 - val_loss: 0.9760
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97601/1 [==============================] - 0s 88ms/step - loss: 0.9760 - val_loss: 0.9727
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97271/1 [==============================] - 0s 79ms/step - loss: 0.9727 - val_loss: 0.9694
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96941/1 [==============================] - 0s 86ms/step - loss: 0.9694 - val_loss: 0.9668
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96681/1 [==============================] - 0s 88ms/step - loss: 0.9668 - val_loss: 0.9645
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96451/1 [==============================] - 0s 81ms/step - loss: 0.9645 - val_loss: 0.9620
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96201/1 [==============================] - 0s 81ms/step - loss: 0.9620 - val_loss: 0.9590
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95901/1 [==============================] - 0s 69ms/step - loss: 0.9590 - val_loss: 0.9563
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95631/1 [==============================] - 0s 69ms/step - loss: 0.9563 - val_loss: 0.9543
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95431/1 [==============================] - 0s 68ms/step - loss: 0.9543 - val_loss: 0.9526
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95261/1 [==============================] - 0s 67ms/step - loss: 0.9526 - val_loss: 0.9506
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95061/1 [==============================] - 0s 69ms/step - loss: 0.9506 - val_loss: 0.9484
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94841/1 [==============================] - 0s 68ms/step - loss: 0.9484 - val_loss: 0.9466
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94661/1 [==============================] - 0s 71ms/step - loss: 0.9466 - val_loss: 0.9448
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94481/1 [==============================] - 0s 69ms/step - loss: 0.9448 - val_loss: 0.9427
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94271/1 [==============================] - 0s 69ms/step - loss: 0.9427 - val_loss: 0.9404
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94041/1 [==============================] - 0s 67ms/step - loss: 0.9404 - val_loss: 0.9383
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93831/1 [==============================] - 0s 83ms/step - loss: 0.9383 - val_loss: 0.9361
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93611/1 [==============================] - 0s 64ms/step - loss: 0.9361 - val_loss: 0.9337
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93371/1 [==============================] - 0s 75ms/step - loss: 0.9337 - val_loss: 0.9312
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93121/1 [==============================] - 0s 72ms/step - loss: 0.9312 - val_loss: 0.9287
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92871/1 [==============================] - 0s 67ms/step - loss: 0.9287 - val_loss: 0.9260
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92601/1 [==============================] - 0s 68ms/step - loss: 0.9260 - val_loss: 0.9231
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92311/1 [==============================] - 0s 70ms/step - loss: 0.9231 - val_loss: 0.9202
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92021/1 [==============================] - 0s 67ms/step - loss: 0.9202 - val_loss: 0.9172
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91721/1 [==============================] - 0s 66ms/step - loss: 0.9172 - val_loss: 0.9140
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91401/1 [==============================] - 0s 69ms/step - loss: 0.9140 - val_loss: 0.9109
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91091/1 [==============================] - 0s 66ms/step - loss: 0.9109 - val_loss: 0.9077
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90771/1 [==============================] - 0s 66ms/step - loss: 0.9077 - val_loss: 0.9044
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90441/1 [==============================] - 0s 81ms/step - loss: 0.9044 - val_loss: 0.9011
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90111/1 [==============================] - 0s 73ms/step - loss: 0.9011 - val_loss: 0.8976
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89761/1 [==============================] - 0s 79ms/step - loss: 0.8976 - val_loss: 0.8941
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89411/1 [==============================] - 0s 81ms/step - loss: 0.8941 - val_loss: 0.8906
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89061/1 [==============================] - 0s 66ms/step - loss: 0.8906 - val_loss: 0.8871
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88711/1 [==============================] - 0s 65ms/step - loss: 0.8871 - val_loss: 0.8835
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88351/1 [==============================] - 0s 72ms/step - loss: 0.8835 - val_loss: 0.8800
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88001/1 [==============================] - 0s 64ms/step - loss: 0.8800 - val_loss: 0.8766
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87661/1 [==============================] - 0s 60ms/step - loss: 0.8766 - val_loss: 0.8731
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87311/1 [==============================] - 0s 59ms/step - loss: 0.8731 - val_loss: 0.8698
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86981/1 [==============================] - 0s 69ms/step - loss: 0.8698 - val_loss: 0.8664
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86641/1 [==============================] - 0s 61ms/step - loss: 0.8664 - val_loss: 0.8632
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86321/1 [==============================] - 0s 82ms/step - loss: 0.8632 - val_loss: 0.8600
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86001/1 [==============================] - 0s 65ms/step - loss: 0.8600 - val_loss: 0.8569
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85691/1 [==============================] - 0s 83ms/step - loss: 0.8569 - val_loss: 0.8538
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85381/1 [==============================] - 0s 70ms/step - loss: 0.8538 - val_loss: 0.8507
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85071/1 [==============================] - 0s 70ms/step - loss: 0.8507 - val_loss: 0.8477
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84771/1 [==============================] - 0s 87ms/step - loss: 0.8477 - val_loss: 0.8448
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84481/1 [==============================] - 0s 80ms/step - loss: 0.8448 - val_loss: 0.8418
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84181/1 [==============================] - 0s 72ms/step - loss: 0.8418 - val_loss: 0.8388
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83881/1 [==============================] - 0s 66ms/step - loss: 0.8388 - val_loss: 0.8359
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83591/1 [==============================] - 0s 85ms/step - loss: 0.8359 - val_loss: 0.8330
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83301/1 [==============================] - 0s 74ms/step - loss: 0.8330 - val_loss: 0.8301
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83011/1 [==============================] - 0s 70ms/step - loss: 0.8301 - val_loss: 0.8273
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82731/1 [==============================] - 0s 72ms/step - loss: 0.8273 - val_loss: 0.8246
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82461/1 [==============================] - 0s 88ms/step - loss: 0.8246 - val_loss: 0.8219
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82191/1 [==============================] - 0s 86ms/step - loss: 0.8219 - val_loss: 0.8190
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81901/1 [==============================] - 0s 83ms/step - loss: 0.8190 - val_loss: 0.8160
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81601/1 [==============================] - 0s 79ms/step - loss: 0.8160 - val_loss: 0.8133
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81331/1 [==============================] - 0s 85ms/step - loss: 0.8133 - val_loss: 0.8108
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81081/1 [==============================] - 0s 93ms/step - loss: 0.8108 - val_loss: 0.8082
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80821/1 [==============================] - 0s 73ms/step - loss: 0.8082 - val_loss: 0.8054
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80541/1 [==============================] - 0s 68ms/step - loss: 0.8054 - val_loss: 0.8027
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80271/1 [==============================] - 0s 69ms/step - loss: 0.8027 - val_loss: 0.8001
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80011/1 [==============================] - 0s 71ms/step - loss: 0.8001 - val_loss: 0.7978
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79781/1 [==============================] - 0s 88ms/step - loss: 0.7978 - val_loss: 0.7953
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79531/1 [==============================] - 0s 75ms/step - loss: 0.7953 - val_loss: 0.7928
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79281/1 [==============================] - 0s 84ms/step - loss: 0.7928 - val_loss: 0.7902
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79021/1 [==============================] - 0s 79ms/step - loss: 0.7902 - val_loss: 0.7878
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78781/1 [==============================] - 0s 70ms/step - loss: 0.7878 - val_loss: 0.7855
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78551/1 [==============================] - 0s 82ms/step - loss: 0.7855 - val_loss: 0.7833
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78331/1 [==============================] - 0s 81ms/step - loss: 0.7833 - val_loss: 0.7810
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78101/1 [==============================] - 0s 81ms/step - loss: 0.7810 - val_loss: 0.7788
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77881/1 [==============================] - 0s 79ms/step - loss: 0.7788 - val_loss: 0.7766
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77661/1 [==============================] - 0s 76ms/step - loss: 0.7766 - val_loss: 0.7744
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77441/1 [==============================] - 0s 72ms/step - loss: 0.7744 - val_loss: 0.7723
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77231/1 [==============================] - 0s 66ms/step - loss: 0.7723 - val_loss: 0.7701
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77011/1 [==============================] - 0s 96ms/step - loss: 0.7701 - val_loss: 0.7679
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76791/1 [==============================] - 0s 91ms/step - loss: 0.7679 - val_loss: 0.7660
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76601/1 [==============================] - 0s 82ms/step - loss: 0.7660 - val_loss: 0.7642
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76421/1 [==============================] - 0s 86ms/step - loss: 0.7642 - val_loss: 0.7623
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76231/1 [==============================] - 0s 95ms/step - loss: 0.7623 - val_loss: 0.7603
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76031/1 [==============================] - 0s 85ms/step - loss: 0.7603 - val_loss: 0.7584
#>  1/61 [..............................] - ETA: 5s44/61 [====================>.........] - ETA: 0s61/61 [==============================] - 0s 1ms/step
#> [1] "[2024-10-23 08:22:37] SE event encoding Finish."
#> [1] "[2024-10-23 08:22:37] step5 Calculate splicing regulation distance and Combine distance ======="
#> [1] "384 rbps are used to calculate splicing regulation information"
#> [1] "Save data"
#> [1] "Save data Finished"
#> [1] "[2024-10-23 08:22:58] step6 Calculate combined event similarity ======="
#> [1] "[2024-10-23 08:22:59] Calculate  A3SS event Similarity"
#> 
#> Attaching package: 'Matrix'
#> The following object is masked from 'package:S4Vectors':
#> 
#>     expand
#> [1] "[2024-10-23 08:23:08] Calculate A3SS event Similarity Finished"
#> [1] "[2024-10-23 08:23:09] Calculate  A5SS event Similarity"
#> [1] "[2024-10-23 08:23:18] Calculate A5SS event Similarity Finished"
#> [1] "[2024-10-23 08:23:19] Calculate  MXE event Similarity"
#> [1] "[2024-10-23 08:23:28] Calculate MXE event Similarity Finished"
#> [1] "[2024-10-23 08:23:29] Calculate  RI event Similarity"
#> [1] "[2024-10-23 08:23:37] Calculate RI event Similarity Finished"
#> [1] "[2024-10-23 08:23:38] Calculate  SE event Similarity"
#> [1] "[2024-10-23 08:23:47] Calculate SE event Similarity Finished"
#> [1] "[2024-10-23 08:23:47] Calculate event similarity Finished."
```

A list of event similarity for different types of splicing events will
be saved to `work_path/imputation/event_similarity/event.similars.rds`;

A list of different types of cell similarity and a list of the number of
neighbors will be saved to rds file
to`work_path/imputation/cell_similarity/cell.similars.rds` and
`work_path/imputation/cell_similarity/dyk.cell.rds`

``` r
print(cellnet.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/cell_similarity/"
cell.similars=readRDS(paste0(cellnet.path,'/cell.similars.rds'))
## three different features can be chosen to quantify cell similarity, including raw event PSI(PSI), and raw junction read counts(RC), and RBP expression(EXP_RBP)
print(names(cell.similars))
#> [1] "PSI"     "RC"      "EXP_RBP"

print(eventnet.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/event_similarity/"
event.similars=readRDS(paste0(eventnet.path,'/event.similars.rds'))
## different splicing event types
print(names(event.similars))
#> [1] "A3SS" "A5SS" "MXE"  "RI"   "SE"
```

#### Parameters used in this step

In this step, some parameters can be adjusted in configure file or the
function parameters directly:

For **cell similarity networks**, SCSES can use RBP expressions, Raw
read count or Raw PSI to measure cell similarities.

`feature_num`: the number of high variable features for PCA before
calculate cell distance.

`rbp`: expression of those RBPs will be used to calculate cell
similarity.

`cell_similarity_data`: data used to calculate cell similarity. Choose
at least on from EXP_RBP, RC, and PSI.

`distance_method`: method used to calculate distance.

`alpha_cell`: restart probability for random walk.

`decay_cell`: threshold of change in the similarity matrix.

`kcell_max`: Maximum number of neighbors.

`kcell_min`: Minimum number of neighbors.

For **event similarity networks**, Event similarities are defined by the
RBP regulatory correlations and an embedding representation by
integrating event sequence similarities.

`ae.para`: parameters of encoding sequence features

`rbp`: expression of those RBPs will be used to calculate RBP regulatory
correlations.

`kevent`: the number of neighbors

`alpha_event`: restart probability for random walk

`decay_event`: threshold of change in the similarity matrix

### Step6. Imputation

Based on these weighted similarity networks, SCSES next will use three
imputation strategies to aggregate the information across similar cells
or events to impute read count or PSI value

``` r
Imputed.data.path = ImputationAll(paras)
#> [1] "[2024-10-23 08:23:48] Get imputed result using cell similarity and event similarity."
#> [1] "Checking data..."
#> [1] "rc checked"
#> [1] "psi checked"
#> [1] "event.info checked"
#> [1] "cell similarity checked"
#> [1] "event similarity checked"
#> [1] "decay=0.05  checked"
#> [1] "Checking event type"
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=PSI;RC;EXP_RBP  checked"
#> Warning in dir.create(output_path):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation' already
#> exists
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/"
#> [1] "[2024-10-23 08:23:49] Running Event_type=A3SS;cell_similarity_feature=PSI"
#> [1] "[2024-10-23 08:23:49] Save data"
#> [1] "[2024-10-23 08:23:49] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500010152.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500010152.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:24:03] Running Event_type=A3SS;cell_similarity_feature=RC"
#> [1] "[2024-10-23 08:24:03] Save data"
#> [1] "[2024-10-23 08:24:03] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-499997995.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-499997995.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:24:16] Running Event_type=A3SS;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-23 08:24:16] Save data"
#> [1] "[2024-10-23 08:24:16] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-499996635.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-499996635.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:24:30] Running Event_type=A5SS;cell_similarity_feature=PSI"
#> [1] "[2024-10-23 08:24:30] Save data"
#> [1] "[2024-10-23 08:24:30] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-499971762.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-499971762.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:24:43] Running Event_type=A5SS;cell_similarity_feature=RC"
#> [1] "[2024-10-23 08:24:43] Save data"
#> [1] "[2024-10-23 08:24:43] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500017092.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500017092.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:24:56] Running Event_type=A5SS;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-23 08:24:56] Save data"
#> [1] "[2024-10-23 08:24:56] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-499980244.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-499980244.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:25:10] Running Event_type=MXE;cell_similarity_feature=PSI"
#> [1] "[2024-10-23 08:25:10] Save data"
#> [1] "[2024-10-23 08:25:10] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500008449.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500008449.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:25:22] Running Event_type=MXE;cell_similarity_feature=RC"
#> [1] "[2024-10-23 08:25:22] Save data"
#> [1] "[2024-10-23 08:25:22] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500003052.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500003052.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:25:35] Running Event_type=MXE;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-23 08:25:35] Save data"
#> [1] "[2024-10-23 08:25:35] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500019338.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500019338.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:25:49] Running Event_type=RI;cell_similarity_feature=PSI"
#> [1] "[2024-10-23 08:25:49] Save data"
#> [1] "[2024-10-23 08:25:49] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500020705.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500020705.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:26:02] Running Event_type=RI;cell_similarity_feature=RC"
#> [1] "[2024-10-23 08:26:02] Save data"
#> [1] "[2024-10-23 08:26:02] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-499971849.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-499971849.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:26:15] Running Event_type=RI;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-23 08:26:15] Save data"
#> [1] "[2024-10-23 08:26:15] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-499977321.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-499977321.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:26:29] Running Event_type=SE;cell_similarity_feature=PSI"
#> [1] "[2024-10-23 08:26:29] Save data"
#> [1] "[2024-10-23 08:26:29] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500003179.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500003179.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:26:42] Running Event_type=SE;cell_similarity_feature=RC"
#> [1] "[2024-10-23 08:26:42] Save data"
#> [1] "[2024-10-23 08:26:42] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500000246.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500000246.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:26:55] Running Event_type=SE;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-23 08:26:55] Save data"
#> [1] "[2024-10-23 08:26:55] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_4082863-500011867.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_4082863-500011867.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-23 08:27:08] Get imputed result using cell similarity and event similarity Finish."
```

Results of each imputation strategy will be saved
into`work_path/imputation`.

``` r
Imputed_seperated = readRDS(Imputed.data.path)
str(Imputed_seperated,max.level=3)
#> List of 2
#>  $ cell      :List of 6
#>   ..$ PSI_PSI    : num [1:2550, 1:15] 0.955 0.829 0.545 0.523 0.284 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ PSI_RC     : num [1:2550, 1:15] 0.889 0.72 0.551 0.82 0.401 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ RC_PSI     : num [1:2550, 1:15] 0.939 0.776 0.635 0.605 0.206 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ RC_RC      : num [1:2550, 1:15] 0.849 0.673 0.563 0.723 0.299 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ EXP_RBP_PSI: num [1:2550, 1:15] 0.908 0.764 0.498 0.582 0.345 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ EXP_RBP_RC : num [1:2550, 1:15] 0.844 0.682 0.513 0.79 0.423 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>  $ cell_event:List of 3
#>   ..$ PSI_PSI    : num [1:2550, 1:15] 0.711 0.633 0.577 0.519 0.595 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ RC_PSI     : num [1:2550, 1:15] 0.676 0.611 0.575 0.509 0.573 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ EXP_RBP_PSI: num [1:2550, 1:15] 0.639 0.603 0.576 0.526 0.578 ...
#>   .. ..- attr(*, "dimnames")=List of 2
# explain:
# For example:
# cell/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw PSI with cell similarities.
# cell/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw inclusion and exclusion read counts with cell similarities, and then calculate the imputed PSI.
# cell_event/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw inclusion and exclusion read counts with cell similarities, and then calculate the imputed PSI. Further impute the results using event similarities.
```

### Step7. Estimation

We recommend different imputation strategies for four scenarios defined
by the abundance of reads counts in the target cell and neighbor cells
(ND, BD, TD+Info, and TD-Info). SCSES pre-trains models to predict the
probability of specific scenario for each cell-event pair. Finally,
SCSES calculates the PSI value using a linear combination of predictions
from the four strategies, weighted by these probabilities.

``` r
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
Imputed.data.final.path = Estimation(paras,rds_imputed_file = Imputed.data.path)
#> [1] "[2024-10-23 08:27:09] Combine imputed psi."
#> [1] "[2024-10-23 08:27:09] Loading data..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_seperated_499982236.rds"
#> [1] "Checking data..."
#> [1] "rc checked"
#> [1] "psi checked"
#> [1] "expr checked"
#> [1] "event checked"
#> [1] "cell similarity checked"
#> [1] "dynamic cell knn checked"
#> [1] "Fine tune model will be used."
#> [1] "classifer checked"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/"
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=PSI;RC;EXP_RBP  checked"
#> [1] "[2024-10-23 08:27:12] Combine imputed psi Finish."
```

A list of final imputation of PSI values will be saved to
`work_path/imputation/Imputed_combined*`.

``` r
Imputed_combined = readRDS(Imputed.data.final.path)
str(Imputed_combined,max.level=2)
#> List of 3
#>  $ PSI    : num [1:2550, 1:15] 0.955 0.829 0.545 0.82 0.284 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ RC     : num [1:2550, 1:15] 0.939 0.776 0.635 0.723 0.206 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ EXP_RBP: num [1:2550, 1:15] 0.908 0.764 0.498 0.79 0.345 ...
#>   ..- attr(*, "dimnames")=List of 2
# The finally imputed PSI were named by cell-cell splicing similarity features, including raw event PSI(PSI), and raw junction read counts(RC), and RBP expression(EXP_RBP).
```

#### Fine-tune the model

To improve the fitness of models for a new dataset, we also provide a
procedure to fine-tune the model. For this analysis, we first build a
reference using a set of splicing events with conserved splicing levels
in different human tissues (<https://zenodo.org/records/6408906>). Then
we compare the splicing level in a new dataset with the reference
records, and give the scenarios definition to each event-cell pair,
which is used to fine-tune the pre-trained model.

The commands to perform these analyses:

``` r
ftrc.path = getFtRawRC(paras)
#> [1] "Loading splicing events for classifer fine tune..."
#> [1] "Checking cells..."
#> [1] "15 cells are considered."
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value_ft/"
#> [1] "[2024-10-23 08:27:12] Counting reads of A3SS events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value_ft//A3SS_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 100%"
#> [1] "[2024-10-23 08:27:15] Counting reads of A3SS events Finish."
#> [1] "[2024-10-23 08:27:15] Counting reads of A5SS events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value_ft//A5SS_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 100%"
#> [1] "[2024-10-23 08:27:17] Counting reads of A5SS events Finish."
#> [1] "[2024-10-23 08:27:17] Counting reads of MXE events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value_ft//MXE_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 100%"
#> [1] "[2024-10-23 08:27:19] Counting reads of MXE events Finish."
#> [1] "[2024-10-23 08:27:19] Counting reads of SE events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value_ft//SE_rjm'
#> already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 100%"
#> [1] "[2024-10-23 08:27:22] Counting reads of SE events Finish."
ftpsi.path = getFtRawPSI(paras)
#> [1] "Checking raw reads..."
#> [1] "Loading splicing events for classifer fine tune..."
#> [1] "[2024-10-23 08:27:22] Calculating PSI value of A3SS events..."
#> [1] "[2024-10-23 08:27:22] Calculating PSI value of A3SS events Finish."
#> [1] "[2024-10-23 08:27:22] Calculating PSI value of A5SS events..."
#> [1] "[2024-10-23 08:27:22] Calculating PSI value of A5SS events Finish."
#> [1] "[2024-10-23 08:27:22] Calculating PSI value of MXE events..."
#> [1] "[2024-10-23 08:27:22] Calculating PSI value of MXE events Finish."
#> [1] "[2024-10-23 08:27:22] Calculating PSI value of SE events..."
#> [1] "[2024-10-23 08:27:22] Calculating PSI value of SE events Finish."
ftrds.path = mergeFtSplicingValue(paras)
ftmodel.path = FtClassifier(paras)
#> [1] "Reading true Ft PSI..."
#> [1] "Loading Pre-training classifer..."
#> [1] "[2024-10-23 08:27:22] Classifer fine tune"
#> [1] "[2024-10-23 08:27:22] Processing raw Ft data..."
#> [1] "Checking data..."
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=PSI;RC;EXP_RBP  checked"
#> [1] "Calculating Classifier Features..."
#> [1] "[2024-10-23 08:27:23] Save data"
#> [1] "[2024-10-23 08:27:23] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-500007830.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-500007830.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:27:46] Save data"
#> [1] "[2024-10-23 08:27:46] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499983124.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499983124.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:28:09] Save data"
#> [1] "[2024-10-23 08:28:09] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-500005175.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-500005175.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:28:33] Save data"
#> [1] "[2024-10-23 08:28:33] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499994832.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499994832.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:28:56] Save data"
#> [1] "[2024-10-23 08:28:56] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499998583.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499998583.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:29:19] Save data"
#> [1] "[2024-10-23 08:29:19] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499977696.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499977696.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:29:42] Save data"
#> [1] "[2024-10-23 08:29:42] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499980757.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499980757.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:30:07] Save data"
#> [1] "[2024-10-23 08:30:07] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-500006628.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-500006628.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:30:32] Save data"
#> [1] "[2024-10-23 08:30:32] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499989448.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499989448.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:30:57] Save data"
#> [1] "[2024-10-23 08:30:57] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499989260.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499989260.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:31:20] Save data"
#> [1] "[2024-10-23 08:31:20] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499968020.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499968020.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:31:44] Save data"
#> [1] "[2024-10-23 08:31:44] Save data Finished"
#> [1] "bash /tmp/RtmpIQFQgd/temp_libpath3e3cee13c6b3e9/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_4082863-499991711.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_4082863-499991711.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-23 08:32:07]  Model training;similarity_type=PSI"
#> [1] "[2024-10-23 08:32:07]  Model training;similarity_type=RC"
#> [1] "[2024-10-23 08:32:08]  Model training;similarity_type=EXP_RBP"
#> [1] "[2024-10-23 08:32:09] Classifer fine tune Finish."
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
ImputedFt.data.final.path = Estimation(paras,rds_imputed_file = Imputed.data.path)
#> [1] "[2024-10-23 08:32:09] Combine imputed psi."
#> [1] "[2024-10-23 08:32:09] Loading data..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_seperated_499982236.rds"
#> [1] "Checking data..."
#> [1] "rc checked"
#> [1] "psi checked"
#> [1] "expr checked"
#> [1] "event checked"
#> [1] "cell similarity checked"
#> [1] "dynamic cell knn checked"
#> [1] "Fine tune model will be used."
#> [1] "classifer checked"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/"
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=PSI;RC;EXP_RBP  checked"
#> [1] "[2024-10-23 08:32:12] Combine imputed psi Finish."

print(ImputedFt.data.final.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_combined_499993925.rds"
```

### Step8. Cell Clustering

``` r
library(umap)
library(ggplot2)
calcu_umap<-function(data,n_neighbors){
  set.seed(12345)
  data = data[which(apply(data, 1, mean) != 0), ]
  data = data[which(apply(data, 1, var) != 0), ]
  D_Reduct_res <- prcomp(t(data), center = T, scale. = T)
  sdev <- D_Reduct_res$sdev
  var_prop <- sdev^2 / sum(sdev^2)
  cumulative_variance <- cumsum(var_prop)
  n_components <- which(cumulative_variance >= 0.7)[1]
  D_Reduct <- D_Reduct_res$x[, 1:n_components]
  umap_res <- umap::umap(D_Reduct,n_neighbors = n_neighbors)
  input <- umap_res$layout
  df<-as.data.frame(input)
  row.names(df)<-colnames(data)
  return(df)
}

mycol=c("#fbb45d","#699ed4","#ef8183")
annotation=read.table("/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/annotation.txt",sep="\t")
Imputed_combined = readRDS(ImputedFt.data.final.path)
data_umap=calcu_umap(Imputed_combined[[3]],n_neighbors = 5)
row.names(data_umap)=gsub(".bam","",row.names(data_umap))
data_umap$group=annotation$V2[match(row.names(data_umap),annotation$V1)]

p=ggplot(data = data_umap,aes(x =V1 ,y =V2))+
    geom_point(aes(fill=group),shape=21,size=1.5,stroke=0.05)+
    scale_fill_manual(values = mycol)+
    xlab("UMAP1")+
    ylab("UMAP2")

print(p)
```

<img src="man/figures/README-unnamed-chunk-18-1.png" width="100%" />

## Example

For iPSC example, see the Rmarkdown tutorials in `analysis/iPSC.Rmd`
