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
#> [1] "[2024-10-22 19:57:07] Detect gene expression: bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/shell/run_featurecounts.sh /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/expr/ /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/test.fa /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/test.gtf /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam 20 cell_line paired /disk/software/subread-2.0.6-source/bin/featureCounts >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/runfeatureCounts.log 2>&1"
#> [1] "[2024-10-22 19:57:34] Detect gene expression Finish."
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
#> [1] "[2024-10-22 19:57:34] Creating Pseudobulk directory..."
#> Warning in dir.create(path = pseudobulk.path, recursive = T):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data' already exists
#> [1] "Pseudobulk bam file exists."
#> [1] "[2024-10-22 19:57:34] Merge Bam Files: /disk/software/samtools/bin/samtools merge -f -@ 20 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data//all.bam /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam/*.bam --no-PG"
#> [1] "[2024-10-22 19:58:16] Merge Bam Files Finish."
#> [1] "[2024-10-22 19:58:16] Bam File Index: /disk/software/samtools/bin/samtools index -@ 20 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data//all.bam"
#> [1] "[2024-10-22 19:58:22] Bam File Index Finish."
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
#> [1] "[2024-10-22 19:58:22] Creating events directory..."
#> Warning in dir.create(path = work_path):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/events' already
#> exists
#> [1] "[2024-10-22 19:58:22] Generating SE event id"
#> arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
#> [1] "[2024-10-22 19:58:22] Generating SE event id Finish."
#> [1] "[2024-10-22 19:58:22] Generating RI event id"
#> arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
#> [1] "[2024-10-22 19:58:37] Generating RI event id Finish."
#> [1] "[2024-10-22 19:58:37] Generating A3SS event id"
#> [1] "[2024-10-22 19:58:37] Generating A3SS event id Finish."
#> [1] "[2024-10-22 19:58:37] Generating A5SS event id"
#> [1] "[2024-10-22 19:58:38] Generating A5SS event id Finish."
#> [1] "[2024-10-22 19:58:38] Generating MXE event id"
#> [1] "[2024-10-22 19:58:38] Generating MXE event id Finish."
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
#> [1] "[2024-10-22 19:58:38] Counting reads of A3SS events..."
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
#> [1] "[2024-10-22 19:58:41] Counting reads of A3SS events Finish."
#> [1] "[2024-10-22 19:58:41] Counting reads of A5SS events..."
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
#> [1] "[2024-10-22 19:58:44] Counting reads of A5SS events Finish."
#> [1] "[2024-10-22 19:58:44] Counting reads of MXE events..."
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
#> [1] "[2024-10-22 19:58:49] Counting reads of MXE events Finish."
#> [1] "[2024-10-22 19:58:49] Counting reads of RI events..."
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
#> [1] "[2024-10-22 19:58:58] Counting reads of RI events Finish."
#> [1] "[2024-10-22 19:58:58] Counting reads of SE events..."
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
#> [1] "[2024-10-22 19:59:28] Counting reads of SE events Finish."
rawpsi.path = getRawPSI(paras)
#> [1] "Checking raw reads..."
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of A3SS events..."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of A3SS events Finish."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of A5SS events..."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of A5SS events Finish."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of MXE events..."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of MXE events Finish."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of RI events..."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of RI events Finish."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of SE events..."
#> [1] "[2024-10-22 19:59:28] Calculating PSI value of SE events Finish."
rawrds.path = mergeSplicingValue(paras)
processed.data.path = preprocessEvent(paras)
#> [1] "[2024-10-22 19:59:29] Processing raw data..."
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
#> [1] "[2024-10-22 19:59:30] Successfully processed data."
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
#> [1] "[2024-10-22 19:59:30] Calculate cell similarity..."
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
#> [1] "[2024-10-22 19:59:30] Computing cell similarity based on PSI"
#> [1] "Calculate similarity among 15 cells."
#> Delta: [0.17081888 0.02261298]
#> [1] "[2024-10-22 20:00:01] Computing cell similarity based on RC"
#> [1] "Calculate similarity among 15 cells."
#> Delta: [0.13512509 0.02119545]
#> [1] "[2024-10-22 20:00:01] Computing cell similarity based on EXP_RBP"
#> [1] "Calculate similarity among 15 cells."
#> [1] "The number of features is greater than the number of rows in the input data."
#> [1] "Total 384 features will be used"
#> Delta: [0.17068148 0.03534559]
#> [1] "[2024-10-22 20:00:01] Calculate cell similarity Finish."
eventnet.path = getEventSimilarity(paras)
#> [1] "[2024-10-22 20:00:01] Calculate events KNN..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds_processed/"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/event_similarity/"
#> [1] "alpha_event=0.8  checked"
#> [1] "kevent=5  checked"
#> [1] "decay_event=0.05  checked"
#> [1] "Checking data..."
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "[2024-10-22 20:00:31] Calculate events feature..."
#> [1] "[2024-10-22 20:00:31] step1 Creating BSgenome for hg19======="
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
#> [1] "[2024-10-22 20:00:31] step2 Extracting features ======="
#> [1] "[2024-10-22 20:00:31] Extracting A3SS features..."
#> [1] "[2024-10-22 20:00:31] Loading events..."
#> [1] "[2024-10-22 20:00:42] Parsing events region..."
#> [1] "[2024-10-22 20:00:51] Extracting length features."
#> [1] "[2024-10-22 20:00:51] Extracting motif features."
#> [1] "[2024-10-22 20:00:51] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-22 20:00:54] Extracting kmer features."
#> [1] "[2024-10-22 20:00:54] Extracting A Ratio features."
#> [1] "[2024-10-22 20:00:56] Saving Result"
#> [1] "[2024-10-22 20:00:56] Extracting A3SS features Finished"
#> [1] "[2024-10-22 20:00:56] Extracting A5SS features..."
#> [1] "[2024-10-22 20:00:56] Loading events..."
#> [1] "[2024-10-22 20:01:07] Parsing events region..."
#> [1] "[2024-10-22 20:01:15] Extracting length features."
#> [1] "[2024-10-22 20:01:15] Extracting motif features."
#> [1] "[2024-10-22 20:01:15] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-22 20:01:17] Extracting kmer features."
#> [1] "[2024-10-22 20:01:17] Extracting A Ratio features."
#> [1] "[2024-10-22 20:01:18] Saving Result"
#> [1] "[2024-10-22 20:01:19] Extracting A5SS features Finished"
#> [1] "[2024-10-22 20:01:19] Extracting MXE features..."
#> [1] "[2024-10-22 20:01:19] Loading events..."
#> [1] "[2024-10-22 20:01:29] Parsing events region..."
#> [1] "[2024-10-22 20:01:37] Extracting length features."
#> [1] "[2024-10-22 20:01:37] Extracting motif features."
#> [1] "[2024-10-22 20:01:37] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-22 20:01:39] Extracting kmer features."
#> [1] "[2024-10-22 20:01:39] Extracting A Ratio features."
#> [1] "[2024-10-22 20:01:39] Saving Result"
#> [1] "[2024-10-22 20:01:39] Extracting MXE features Finished"
#> [1] "[2024-10-22 20:01:39] Extracting RI features..."
#> [1] "[2024-10-22 20:01:39] Loading events..."
#> [1] "[2024-10-22 20:01:49] Parsing events region..."
#> [1] "[2024-10-22 20:01:58] Extracting length features."
#> [1] "[2024-10-22 20:01:58] Extracting motif features."
#> [1] "[2024-10-22 20:01:58] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-22 20:02:01] Extracting kmer features."
#> [1] "[2024-10-22 20:02:01] Extracting A Ratio features."
#> [1] "[2024-10-22 20:02:02] Saving Result"
#> [1] "[2024-10-22 20:02:02] Extracting RI features Finished"
#> [1] "[2024-10-22 20:02:02] Extracting SE features..."
#> [1] "[2024-10-22 20:02:02] Loading events..."
#> [1] "[2024-10-22 20:02:13] Parsing events region..."
#> [1] "[2024-10-22 20:02:26] Extracting length features."
#> [1] "[2024-10-22 20:02:29] Extracting motif features."
#> [1] "[2024-10-22 20:02:29] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-10-22 20:02:41] Extracting kmer features."
#> [1] "[2024-10-22 20:02:41] Extracting A Ratio features."
#> [1] "[2024-10-22 20:02:41] Saving Result"
#> [1] "[2024-10-22 20:02:42] Extracting SE features Finished"
#> [1] "[2024-10-22 20:02:42] step3 Combining events feature ======="
#> [1] "[2024-10-22 20:02:42] Parsing A3SS features..."
#> [1] "[2024-10-22 20:02:42] Parsing A3SS features Finished"
#> [1] "[2024-10-22 20:02:42] Parsing A5SS features..."
#> [1] "[2024-10-22 20:02:42] Parsing A5SS features Finished"
#> [1] "[2024-10-22 20:02:42] Parsing MXE features..."
#> [1] "[2024-10-22 20:02:42] Parsing MXE features Finished"
#> [1] "[2024-10-22 20:02:42] Parsing RI features..."
#> [1] "[2024-10-22 20:02:43] Parsing RI features Finished"
#> [1] "[2024-10-22 20:02:43] Parsing SE features..."
#> [1] "[2024-10-22 20:02:45] Parsing SE features Finished"
#> [1] "[2024-10-22 20:02:45] step4 Encoding events feature ======="
#> [1] "[2024-10-22 20:02:45] A3SS event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02241/1 [==============================] - 2s 2s/step - loss: 1.0224 - val_loss: 1.0158
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01581/1 [==============================] - 0s 78ms/step - loss: 1.0158 - val_loss: 1.0133
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01331/1 [==============================] - 0s 67ms/step - loss: 1.0133 - val_loss: 1.0112
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01121/1 [==============================] - 0s 57ms/step - loss: 1.0112 - val_loss: 1.0090
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00901/1 [==============================] - 0s 56ms/step - loss: 1.0090 - val_loss: 1.0061
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00611/1 [==============================] - 0s 54ms/step - loss: 1.0061 - val_loss: 1.0021
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00211/1 [==============================] - 0s 71ms/step - loss: 1.0021 - val_loss: 0.9966
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99661/1 [==============================] - 0s 56ms/step - loss: 0.9966 - val_loss: 0.9896
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98961/1 [==============================] - 0s 69ms/step - loss: 0.9896 - val_loss: 0.9817
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98171/1 [==============================] - 0s 71ms/step - loss: 0.9817 - val_loss: 0.9749
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97491/1 [==============================] - 0s 74ms/step - loss: 0.9749 - val_loss: 0.9674
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96741/1 [==============================] - 0s 72ms/step - loss: 0.9674 - val_loss: 0.9585
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95851/1 [==============================] - 0s 79ms/step - loss: 0.9585 - val_loss: 0.9504
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95041/1 [==============================] - 0s 56ms/step - loss: 0.9504 - val_loss: 0.9432
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94321/1 [==============================] - 0s 50ms/step - loss: 0.9432 - val_loss: 0.9360
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93601/1 [==============================] - 0s 51ms/step - loss: 0.9360 - val_loss: 0.9290
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92901/1 [==============================] - 0s 49ms/step - loss: 0.9290 - val_loss: 0.9221
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92211/1 [==============================] - 0s 51ms/step - loss: 0.9221 - val_loss: 0.9145
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91451/1 [==============================] - 0s 69ms/step - loss: 0.9145 - val_loss: 0.9065
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90651/1 [==============================] - 0s 54ms/step - loss: 0.9065 - val_loss: 0.8981
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89811/1 [==============================] - 0s 72ms/step - loss: 0.8981 - val_loss: 0.8890
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88901/1 [==============================] - 0s 56ms/step - loss: 0.8890 - val_loss: 0.8793
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87931/1 [==============================] - 0s 69ms/step - loss: 0.8793 - val_loss: 0.8690
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86901/1 [==============================] - 0s 71ms/step - loss: 0.8690 - val_loss: 0.8584
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85841/1 [==============================] - 0s 73ms/step - loss: 0.8584 - val_loss: 0.8477
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84771/1 [==============================] - 0s 69ms/step - loss: 0.8477 - val_loss: 0.8369
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83691/1 [==============================] - 0s 54ms/step - loss: 0.8369 - val_loss: 0.8259
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82591/1 [==============================] - 0s 54ms/step - loss: 0.8259 - val_loss: 0.8148
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81481/1 [==============================] - 0s 70ms/step - loss: 0.8148 - val_loss: 0.8037
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80371/1 [==============================] - 0s 54ms/step - loss: 0.8037 - val_loss: 0.7927
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79271/1 [==============================] - 0s 50ms/step - loss: 0.7927 - val_loss: 0.7818
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78181/1 [==============================] - 0s 59ms/step - loss: 0.7818 - val_loss: 0.7709
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77091/1 [==============================] - 0s 61ms/step - loss: 0.7709 - val_loss: 0.7602
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76021/1 [==============================] - 0s 63ms/step - loss: 0.7602 - val_loss: 0.7496
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74961/1 [==============================] - 0s 54ms/step - loss: 0.7496 - val_loss: 0.7391
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73911/1 [==============================] - 0s 52ms/step - loss: 0.7391 - val_loss: 0.7286
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72861/1 [==============================] - 0s 51ms/step - loss: 0.7286 - val_loss: 0.7183
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71831/1 [==============================] - 0s 58ms/step - loss: 0.7183 - val_loss: 0.7081
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70811/1 [==============================] - 0s 52ms/step - loss: 0.7081 - val_loss: 0.6981
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69811/1 [==============================] - 0s 53ms/step - loss: 0.6981 - val_loss: 0.6882
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68821/1 [==============================] - 0s 72ms/step - loss: 0.6882 - val_loss: 0.6785
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67851/1 [==============================] - 0s 55ms/step - loss: 0.6785 - val_loss: 0.6688
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66881/1 [==============================] - 0s 84ms/step - loss: 0.6688 - val_loss: 0.6592
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65921/1 [==============================] - 0s 77ms/step - loss: 0.6592 - val_loss: 0.6497
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64971/1 [==============================] - 0s 59ms/step - loss: 0.6497 - val_loss: 0.6402
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64021/1 [==============================] - 0s 74ms/step - loss: 0.6402 - val_loss: 0.6308
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63081/1 [==============================] - 0s 56ms/step - loss: 0.6308 - val_loss: 0.6214
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62141/1 [==============================] - 0s 54ms/step - loss: 0.6214 - val_loss: 0.6122
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61221/1 [==============================] - 0s 53ms/step - loss: 0.6122 - val_loss: 0.6031
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60311/1 [==============================] - 0s 56ms/step - loss: 0.6031 - val_loss: 0.5940
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.59401/1 [==============================] - 0s 68ms/step - loss: 0.5940 - val_loss: 0.5851
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.58511/1 [==============================] - 0s 56ms/step - loss: 0.5851 - val_loss: 0.5763
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.57631/1 [==============================] - 0s 51ms/step - loss: 0.5763 - val_loss: 0.5675
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.56751/1 [==============================] - 0s 51ms/step - loss: 0.5675 - val_loss: 0.5588
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.55881/1 [==============================] - 0s 54ms/step - loss: 0.5588 - val_loss: 0.5503
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.55031/1 [==============================] - 0s 72ms/step - loss: 0.5503 - val_loss: 0.5418
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.54181/1 [==============================] - 0s 56ms/step - loss: 0.5418 - val_loss: 0.5334
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.53341/1 [==============================] - 0s 57ms/step - loss: 0.5334 - val_loss: 0.5251
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.52511/1 [==============================] - 0s 74ms/step - loss: 0.5251 - val_loss: 0.5169
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.51691/1 [==============================] - 0s 59ms/step - loss: 0.5169 - val_loss: 0.5088
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.50881/1 [==============================] - 0s 87ms/step - loss: 0.5088 - val_loss: 0.5007
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.50071/1 [==============================] - 0s 70ms/step - loss: 0.5007 - val_loss: 0.4927
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.49271/1 [==============================] - 0s 72ms/step - loss: 0.4927 - val_loss: 0.4849
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.48491/1 [==============================] - 0s 67ms/step - loss: 0.4849 - val_loss: 0.4771
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.47711/1 [==============================] - 0s 78ms/step - loss: 0.4771 - val_loss: 0.4694
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.46941/1 [==============================] - 0s 96ms/step - loss: 0.4694 - val_loss: 0.4617
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.46171/1 [==============================] - 0s 78ms/step - loss: 0.4617 - val_loss: 0.4542
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.45421/1 [==============================] - 0s 61ms/step - loss: 0.4542 - val_loss: 0.4467
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.44671/1 [==============================] - 0s 67ms/step - loss: 0.4467 - val_loss: 0.4393
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.43931/1 [==============================] - 0s 90ms/step - loss: 0.4393 - val_loss: 0.4319
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.43191/1 [==============================] - 0s 90ms/step - loss: 0.4319 - val_loss: 0.4247
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.42471/1 [==============================] - 0s 76ms/step - loss: 0.4247 - val_loss: 0.4175
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.41751/1 [==============================] - 0s 57ms/step - loss: 0.4175 - val_loss: 0.4104
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.41041/1 [==============================] - 0s 55ms/step - loss: 0.4104 - val_loss: 0.4034
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.40341/1 [==============================] - 0s 62ms/step - loss: 0.4034 - val_loss: 0.3964
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.39641/1 [==============================] - 0s 71ms/step - loss: 0.3964 - val_loss: 0.3895
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.38951/1 [==============================] - 0s 58ms/step - loss: 0.3895 - val_loss: 0.3828
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.38281/1 [==============================] - 0s 75ms/step - loss: 0.3828 - val_loss: 0.3760
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.37601/1 [==============================] - 0s 71ms/step - loss: 0.3760 - val_loss: 0.3694
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36941/1 [==============================] - 0s 78ms/step - loss: 0.3694 - val_loss: 0.3629
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36291/1 [==============================] - 0s 77ms/step - loss: 0.3629 - val_loss: 0.3564
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.35641/1 [==============================] - 0s 62ms/step - loss: 0.3564 - val_loss: 0.3500
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.35001/1 [==============================] - 0s 69ms/step - loss: 0.3500 - val_loss: 0.3437
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.34371/1 [==============================] - 0s 74ms/step - loss: 0.3437 - val_loss: 0.3376
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.33761/1 [==============================] - 0s 55ms/step - loss: 0.3376 - val_loss: 0.3314
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.33141/1 [==============================] - 0s 76ms/step - loss: 0.3314 - val_loss: 0.3254
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.32541/1 [==============================] - 0s 76ms/step - loss: 0.3254 - val_loss: 0.3195
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.31951/1 [==============================] - 0s 77ms/step - loss: 0.3195 - val_loss: 0.3137
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.31371/1 [==============================] - 0s 68ms/step - loss: 0.3137 - val_loss: 0.3079
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30791/1 [==============================] - 0s 73ms/step - loss: 0.3079 - val_loss: 0.3024
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30241/1 [==============================] - 0s 78ms/step - loss: 0.3024 - val_loss: 0.2968
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29681/1 [==============================] - 0s 71ms/step - loss: 0.2968 - val_loss: 0.2914
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29141/1 [==============================] - 0s 82ms/step - loss: 0.2914 - val_loss: 0.2861
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.28611/1 [==============================] - 0s 72ms/step - loss: 0.2861 - val_loss: 0.2808
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.28081/1 [==============================] - 0s 59ms/step - loss: 0.2808 - val_loss: 0.2756
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27561/1 [==============================] - 0s 74ms/step - loss: 0.2756 - val_loss: 0.2705
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27051/1 [==============================] - 0s 65ms/step - loss: 0.2705 - val_loss: 0.2656
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.26561/1 [==============================] - 0s 77ms/step - loss: 0.2656 - val_loss: 0.2607
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.26071/1 [==============================] - 0s 72ms/step - loss: 0.2607 - val_loss: 0.2559
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25591/1 [==============================] - 0s 56ms/step - loss: 0.2559 - val_loss: 0.2512
#> 1/4 [======>.......................] - ETA: 0s4/4 [==============================] - 0s 1ms/step
#> [1] "[2024-10-22 20:02:54] A3SS event encoding Finish."
#> [1] "[2024-10-22 20:02:54] A5SS event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02601/1 [==============================] - 1s 1s/step - loss: 1.0260 - val_loss: 1.0165
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01651/1 [==============================] - 0s 57ms/step - loss: 1.0165 - val_loss: 1.0129
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01291/1 [==============================] - 0s 78ms/step - loss: 1.0129 - val_loss: 1.0108
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01081/1 [==============================] - 0s 63ms/step - loss: 1.0108 - val_loss: 1.0090
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00901/1 [==============================] - 0s 66ms/step - loss: 1.0090 - val_loss: 1.0071
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00711/1 [==============================] - 0s 74ms/step - loss: 1.0071 - val_loss: 1.0049
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00491/1 [==============================] - 0s 61ms/step - loss: 1.0049 - val_loss: 1.0020
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00201/1 [==============================] - 0s 72ms/step - loss: 1.0020 - val_loss: 0.9981
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99811/1 [==============================] - 0s 55ms/step - loss: 0.9981 - val_loss: 0.9928
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99281/1 [==============================] - 0s 72ms/step - loss: 0.9928 - val_loss: 0.9860
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98601/1 [==============================] - 0s 59ms/step - loss: 0.9860 - val_loss: 0.9776
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97761/1 [==============================] - 0s 66ms/step - loss: 0.9776 - val_loss: 0.9687
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96871/1 [==============================] - 0s 73ms/step - loss: 0.9687 - val_loss: 0.9602
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96021/1 [==============================] - 0s 81ms/step - loss: 0.9602 - val_loss: 0.9501
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95011/1 [==============================] - 0s 74ms/step - loss: 0.9501 - val_loss: 0.9380
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93801/1 [==============================] - 0s 71ms/step - loss: 0.9380 - val_loss: 0.9261
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92611/1 [==============================] - 0s 56ms/step - loss: 0.9261 - val_loss: 0.9146
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91461/1 [==============================] - 0s 92ms/step - loss: 0.9146 - val_loss: 0.9026
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90261/1 [==============================] - 0s 57ms/step - loss: 0.9026 - val_loss: 0.8908
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89081/1 [==============================] - 0s 72ms/step - loss: 0.8908 - val_loss: 0.8797
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87971/1 [==============================] - 0s 68ms/step - loss: 0.8797 - val_loss: 0.8681
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86811/1 [==============================] - 0s 70ms/step - loss: 0.8681 - val_loss: 0.8561
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85611/1 [==============================] - 0s 75ms/step - loss: 0.8561 - val_loss: 0.8445
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84451/1 [==============================] - 0s 71ms/step - loss: 0.8445 - val_loss: 0.8332
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83321/1 [==============================] - 0s 56ms/step - loss: 0.8332 - val_loss: 0.8219
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82191/1 [==============================] - 0s 52ms/step - loss: 0.8219 - val_loss: 0.8108
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81081/1 [==============================] - 0s 76ms/step - loss: 0.8108 - val_loss: 0.7999
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79991/1 [==============================] - 0s 57ms/step - loss: 0.7999 - val_loss: 0.7889
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78891/1 [==============================] - 0s 54ms/step - loss: 0.7889 - val_loss: 0.7779
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77791/1 [==============================] - 0s 52ms/step - loss: 0.7779 - val_loss: 0.7669
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76691/1 [==============================] - 0s 70ms/step - loss: 0.7669 - val_loss: 0.7559
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75591/1 [==============================] - 0s 57ms/step - loss: 0.7559 - val_loss: 0.7449
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74491/1 [==============================] - 0s 66ms/step - loss: 0.7449 - val_loss: 0.7338
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73381/1 [==============================] - 0s 56ms/step - loss: 0.7338 - val_loss: 0.7228
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72281/1 [==============================] - 0s 53ms/step - loss: 0.7228 - val_loss: 0.7119
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71191/1 [==============================] - 0s 72ms/step - loss: 0.7119 - val_loss: 0.7009
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70091/1 [==============================] - 0s 55ms/step - loss: 0.7009 - val_loss: 0.6900
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69001/1 [==============================] - 0s 55ms/step - loss: 0.6900 - val_loss: 0.6790
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67901/1 [==============================] - 0s 68ms/step - loss: 0.6790 - val_loss: 0.6682
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66821/1 [==============================] - 0s 53ms/step - loss: 0.6682 - val_loss: 0.6574
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65741/1 [==============================] - 0s 80ms/step - loss: 0.6574 - val_loss: 0.6467
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64671/1 [==============================] - 0s 61ms/step - loss: 0.6467 - val_loss: 0.6360
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63601/1 [==============================] - 0s 70ms/step - loss: 0.6360 - val_loss: 0.6255
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62551/1 [==============================] - 0s 55ms/step - loss: 0.6255 - val_loss: 0.6151
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61511/1 [==============================] - 0s 58ms/step - loss: 0.6151 - val_loss: 0.6048
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60481/1 [==============================] - 0s 75ms/step - loss: 0.6048 - val_loss: 0.5946
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.59461/1 [==============================] - 0s 62ms/step - loss: 0.5946 - val_loss: 0.5845
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.58451/1 [==============================] - 0s 75ms/step - loss: 0.5845 - val_loss: 0.5745
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.57451/1 [==============================] - 0s 71ms/step - loss: 0.5745 - val_loss: 0.5646
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.56461/1 [==============================] - 0s 75ms/step - loss: 0.5646 - val_loss: 0.5549
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.55491/1 [==============================] - 0s 73ms/step - loss: 0.5549 - val_loss: 0.5452
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.54521/1 [==============================] - 0s 68ms/step - loss: 0.5452 - val_loss: 0.5356
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.53561/1 [==============================] - 0s 68ms/step - loss: 0.5356 - val_loss: 0.5262
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.52621/1 [==============================] - 0s 57ms/step - loss: 0.5262 - val_loss: 0.5169
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.51691/1 [==============================] - 0s 73ms/step - loss: 0.5169 - val_loss: 0.5076
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.50761/1 [==============================] - 0s 66ms/step - loss: 0.5076 - val_loss: 0.4985
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.49851/1 [==============================] - 0s 78ms/step - loss: 0.4985 - val_loss: 0.4895
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.48951/1 [==============================] - 0s 55ms/step - loss: 0.4895 - val_loss: 0.4806
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.48061/1 [==============================] - 0s 54ms/step - loss: 0.4806 - val_loss: 0.4718
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.47181/1 [==============================] - 0s 67ms/step - loss: 0.4718 - val_loss: 0.4632
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.46321/1 [==============================] - 0s 76ms/step - loss: 0.4632 - val_loss: 0.4547
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.45471/1 [==============================] - 0s 73ms/step - loss: 0.4547 - val_loss: 0.4463
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.44631/1 [==============================] - 0s 56ms/step - loss: 0.4463 - val_loss: 0.4380
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.43801/1 [==============================] - 0s 77ms/step - loss: 0.4380 - val_loss: 0.4298
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.42981/1 [==============================] - 0s 63ms/step - loss: 0.4298 - val_loss: 0.4218
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.42181/1 [==============================] - 0s 54ms/step - loss: 0.4218 - val_loss: 0.4139
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.41391/1 [==============================] - 0s 83ms/step - loss: 0.4139 - val_loss: 0.4061
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.40611/1 [==============================] - 0s 72ms/step - loss: 0.4061 - val_loss: 0.3984
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.39841/1 [==============================] - 0s 67ms/step - loss: 0.3984 - val_loss: 0.3908
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.39081/1 [==============================] - 0s 55ms/step - loss: 0.3908 - val_loss: 0.3833
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.38331/1 [==============================] - 0s 66ms/step - loss: 0.3833 - val_loss: 0.3760
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.37601/1 [==============================] - 0s 76ms/step - loss: 0.3760 - val_loss: 0.3688
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36881/1 [==============================] - 0s 66ms/step - loss: 0.3688 - val_loss: 0.3617
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36171/1 [==============================] - 0s 72ms/step - loss: 0.3617 - val_loss: 0.3546
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.35461/1 [==============================] - 0s 57ms/step - loss: 0.3546 - val_loss: 0.3477
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.34771/1 [==============================] - 0s 75ms/step - loss: 0.3477 - val_loss: 0.3409
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.34091/1 [==============================] - 0s 72ms/step - loss: 0.3409 - val_loss: 0.3342
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.33421/1 [==============================] - 0s 59ms/step - loss: 0.3342 - val_loss: 0.3276
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.32761/1 [==============================] - 0s 69ms/step - loss: 0.3276 - val_loss: 0.3211
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.32111/1 [==============================] - 0s 74ms/step - loss: 0.3211 - val_loss: 0.3148
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.31481/1 [==============================] - 0s 70ms/step - loss: 0.3148 - val_loss: 0.3085
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30851/1 [==============================] - 0s 69ms/step - loss: 0.3085 - val_loss: 0.3023
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30231/1 [==============================] - 0s 55ms/step - loss: 0.3023 - val_loss: 0.2962
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29621/1 [==============================] - 0s 53ms/step - loss: 0.2962 - val_loss: 0.2902
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29021/1 [==============================] - 0s 54ms/step - loss: 0.2902 - val_loss: 0.2842
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.28421/1 [==============================] - 0s 74ms/step - loss: 0.2842 - val_loss: 0.2784
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27841/1 [==============================] - 0s 62ms/step - loss: 0.2784 - val_loss: 0.2727
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27271/1 [==============================] - 0s 73ms/step - loss: 0.2727 - val_loss: 0.2670
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.26701/1 [==============================] - 0s 57ms/step - loss: 0.2670 - val_loss: 0.2615
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.26151/1 [==============================] - 0s 51ms/step - loss: 0.2615 - val_loss: 0.2560
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25601/1 [==============================] - 0s 52ms/step - loss: 0.2560 - val_loss: 0.2507
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25071/1 [==============================] - 0s 75ms/step - loss: 0.2507 - val_loss: 0.2454
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.24541/1 [==============================] - 0s 71ms/step - loss: 0.2454 - val_loss: 0.2403
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.24031/1 [==============================] - 0s 68ms/step - loss: 0.2403 - val_loss: 0.2352
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.23521/1 [==============================] - 0s 71ms/step - loss: 0.2352 - val_loss: 0.2303
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.23031/1 [==============================] - 0s 60ms/step - loss: 0.2303 - val_loss: 0.2254
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.22541/1 [==============================] - 0s 64ms/step - loss: 0.2254 - val_loss: 0.2206
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.22061/1 [==============================] - 0s 64ms/step - loss: 0.2206 - val_loss: 0.2159
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.21591/1 [==============================] - 0s 67ms/step - loss: 0.2159 - val_loss: 0.2113
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.21131/1 [==============================] - 0s 54ms/step - loss: 0.2113 - val_loss: 0.2069
#> 1/4 [======>.......................] - ETA: 0s4/4 [==============================] - 0s 1ms/step
#> [1] "[2024-10-22 20:03:03] A5SS event encoding Finish."
#> [1] "[2024-10-22 20:03:03] MXE event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98541/1 [==============================] - 1s 1s/step - loss: 0.9854 - val_loss: 0.9763
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97631/1 [==============================] - 0s 51ms/step - loss: 0.9763 - val_loss: 0.9728
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97281/1 [==============================] - 0s 66ms/step - loss: 0.9728 - val_loss: 0.9702
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97021/1 [==============================] - 0s 50ms/step - loss: 0.9702 - val_loss: 0.9674
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96741/1 [==============================] - 0s 68ms/step - loss: 0.9674 - val_loss: 0.9638
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96381/1 [==============================] - 0s 66ms/step - loss: 0.9638 - val_loss: 0.9592
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95921/1 [==============================] - 0s 68ms/step - loss: 0.9592 - val_loss: 0.9532
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95321/1 [==============================] - 0s 70ms/step - loss: 0.9532 - val_loss: 0.9452
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94521/1 [==============================] - 0s 51ms/step - loss: 0.9452 - val_loss: 0.9353
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93531/1 [==============================] - 0s 47ms/step - loss: 0.9353 - val_loss: 0.9238
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92381/1 [==============================] - 0s 67ms/step - loss: 0.9238 - val_loss: 0.9116
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91161/1 [==============================] - 0s 60ms/step - loss: 0.9116 - val_loss: 0.8979
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89791/1 [==============================] - 0s 68ms/step - loss: 0.8979 - val_loss: 0.8822
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88221/1 [==============================] - 0s 53ms/step - loss: 0.8822 - val_loss: 0.8656
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86561/1 [==============================] - 0s 45ms/step - loss: 0.8656 - val_loss: 0.8486
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84861/1 [==============================] - 0s 46ms/step - loss: 0.8486 - val_loss: 0.8308
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83081/1 [==============================] - 0s 46ms/step - loss: 0.8308 - val_loss: 0.8120
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81201/1 [==============================] - 0s 46ms/step - loss: 0.8120 - val_loss: 0.7930
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79301/1 [==============================] - 0s 51ms/step - loss: 0.7930 - val_loss: 0.7736
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77361/1 [==============================] - 0s 48ms/step - loss: 0.7736 - val_loss: 0.7540
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75401/1 [==============================] - 0s 65ms/step - loss: 0.7540 - val_loss: 0.7344
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73441/1 [==============================] - 0s 57ms/step - loss: 0.7344 - val_loss: 0.7148
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71481/1 [==============================] - 0s 63ms/step - loss: 0.7148 - val_loss: 0.6953
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69531/1 [==============================] - 0s 60ms/step - loss: 0.6953 - val_loss: 0.6761
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67611/1 [==============================] - 0s 65ms/step - loss: 0.6761 - val_loss: 0.6571
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65711/1 [==============================] - 0s 51ms/step - loss: 0.6571 - val_loss: 0.6383
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63831/1 [==============================] - 0s 48ms/step - loss: 0.6383 - val_loss: 0.6198
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61981/1 [==============================] - 0s 60ms/step - loss: 0.6198 - val_loss: 0.6016
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60161/1 [==============================] - 0s 58ms/step - loss: 0.6016 - val_loss: 0.5837
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.58371/1 [==============================] - 0s 49ms/step - loss: 0.5837 - val_loss: 0.5660
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.56601/1 [==============================] - 0s 62ms/step - loss: 0.5660 - val_loss: 0.5486
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.54861/1 [==============================] - 0s 52ms/step - loss: 0.5486 - val_loss: 0.5315
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.53151/1 [==============================] - 0s 47ms/step - loss: 0.5315 - val_loss: 0.5149
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.51491/1 [==============================] - 0s 47ms/step - loss: 0.5149 - val_loss: 0.4986
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.49861/1 [==============================] - 0s 51ms/step - loss: 0.4986 - val_loss: 0.4826
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.48261/1 [==============================] - 0s 71ms/step - loss: 0.4826 - val_loss: 0.4670
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.46701/1 [==============================] - 0s 61ms/step - loss: 0.4670 - val_loss: 0.4518
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.45181/1 [==============================] - 0s 57ms/step - loss: 0.4518 - val_loss: 0.4369
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.43691/1 [==============================] - 0s 57ms/step - loss: 0.4369 - val_loss: 0.4223
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.42231/1 [==============================] - 0s 51ms/step - loss: 0.4223 - val_loss: 0.4080
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.40801/1 [==============================] - 0s 50ms/step - loss: 0.4080 - val_loss: 0.3939
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.39391/1 [==============================] - 0s 52ms/step - loss: 0.3939 - val_loss: 0.3802
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.38021/1 [==============================] - 0s 52ms/step - loss: 0.3802 - val_loss: 0.3667
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.36671/1 [==============================] - 0s 64ms/step - loss: 0.3667 - val_loss: 0.3534
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.35341/1 [==============================] - 0s 53ms/step - loss: 0.3534 - val_loss: 0.3405
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.34051/1 [==============================] - 0s 72ms/step - loss: 0.3405 - val_loss: 0.3277
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.32771/1 [==============================] - 0s 69ms/step - loss: 0.3277 - val_loss: 0.3153
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.31531/1 [==============================] - 0s 65ms/step - loss: 0.3153 - val_loss: 0.3030
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.30301/1 [==============================] - 0s 64ms/step - loss: 0.3030 - val_loss: 0.2911
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.29111/1 [==============================] - 0s 70ms/step - loss: 0.2911 - val_loss: 0.2795
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.27951/1 [==============================] - 0s 65ms/step - loss: 0.2795 - val_loss: 0.2681
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.26811/1 [==============================] - 0s 55ms/step - loss: 0.2681 - val_loss: 0.2571
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.25711/1 [==============================] - 0s 71ms/step - loss: 0.2571 - val_loss: 0.2463
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.24631/1 [==============================] - 0s 56ms/step - loss: 0.2463 - val_loss: 0.2358
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.23581/1 [==============================] - 0s 76ms/step - loss: 0.2358 - val_loss: 0.2257
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.22571/1 [==============================] - 0s 71ms/step - loss: 0.2257 - val_loss: 0.2159
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.21591/1 [==============================] - 0s 69ms/step - loss: 0.2159 - val_loss: 0.2063
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.20631/1 [==============================] - 0s 61ms/step - loss: 0.2063 - val_loss: 0.1971
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.19711/1 [==============================] - 0s 71ms/step - loss: 0.1971 - val_loss: 0.1882
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.18821/1 [==============================] - 0s 70ms/step - loss: 0.1882 - val_loss: 0.1797
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.17971/1 [==============================] - 0s 55ms/step - loss: 0.1797 - val_loss: 0.1715
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.17151/1 [==============================] - 0s 68ms/step - loss: 0.1715 - val_loss: 0.1636
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.16361/1 [==============================] - 0s 53ms/step - loss: 0.1636 - val_loss: 0.1560
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.15601/1 [==============================] - 0s 57ms/step - loss: 0.1560 - val_loss: 0.1488
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.14881/1 [==============================] - 0s 63ms/step - loss: 0.1488 - val_loss: 0.1419
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.14191/1 [==============================] - 0s 65ms/step - loss: 0.1419 - val_loss: 0.1353
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.13531/1 [==============================] - 0s 69ms/step - loss: 0.1353 - val_loss: 0.1290
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.12901/1 [==============================] - 0s 69ms/step - loss: 0.1290 - val_loss: 0.1231
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.12311/1 [==============================] - 0s 66ms/step - loss: 0.1231 - val_loss: 0.1175
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.11751/1 [==============================] - 0s 68ms/step - loss: 0.1175 - val_loss: 0.1122
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.11221/1 [==============================] - 0s 52ms/step - loss: 0.1122 - val_loss: 0.1073
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.10731/1 [==============================] - 0s 47ms/step - loss: 0.1073 - val_loss: 0.1026
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.10261/1 [==============================] - 0s 59ms/step - loss: 0.1026 - val_loss: 0.0983
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.09831/1 [==============================] - 0s 58ms/step - loss: 0.0983 - val_loss: 0.0942
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.09421/1 [==============================] - 0s 74ms/step - loss: 0.0942 - val_loss: 0.0904
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.09041/1 [==============================] - 0s 64ms/step - loss: 0.0904 - val_loss: 0.0868
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.08681/1 [==============================] - 0s 66ms/step - loss: 0.0868 - val_loss: 0.0834
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.08341/1 [==============================] - 0s 68ms/step - loss: 0.0834 - val_loss: 0.0802
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.08021/1 [==============================] - 0s 54ms/step - loss: 0.0802 - val_loss: 0.0773
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.07731/1 [==============================] - 0s 47ms/step - loss: 0.0773 - val_loss: 0.0745
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.07451/1 [==============================] - 0s 68ms/step - loss: 0.0745 - val_loss: 0.0718
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.07181/1 [==============================] - 0s 56ms/step - loss: 0.0718 - val_loss: 0.0693
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06931/1 [==============================] - 0s 51ms/step - loss: 0.0693 - val_loss: 0.0670
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06701/1 [==============================] - 0s 69ms/step - loss: 0.0670 - val_loss: 0.0648
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06481/1 [==============================] - 0s 59ms/step - loss: 0.0648 - val_loss: 0.0628
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06281/1 [==============================] - 0s 49ms/step - loss: 0.0628 - val_loss: 0.0609
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.06091/1 [==============================] - 0s 52ms/step - loss: 0.0609 - val_loss: 0.0591
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05911/1 [==============================] - 0s 72ms/step - loss: 0.0591 - val_loss: 0.0575
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05751/1 [==============================] - 0s 72ms/step - loss: 0.0575 - val_loss: 0.0559
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05591/1 [==============================] - 0s 68ms/step - loss: 0.0559 - val_loss: 0.0545
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05451/1 [==============================] - 0s 60ms/step - loss: 0.0545 - val_loss: 0.0532
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05321/1 [==============================] - 0s 75ms/step - loss: 0.0532 - val_loss: 0.0520
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05201/1 [==============================] - 0s 66ms/step - loss: 0.0520 - val_loss: 0.0508
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.05081/1 [==============================] - 0s 68ms/step - loss: 0.0508 - val_loss: 0.0497
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.04971/1 [==============================] - 0s 67ms/step - loss: 0.0497 - val_loss: 0.0486
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.04861/1 [==============================] - 0s 69ms/step - loss: 0.0486 - val_loss: 0.0477
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.04771/1 [==============================] - 0s 68ms/step - loss: 0.0477 - val_loss: 0.0468
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.04681/1 [==============================] - 0s 51ms/step - loss: 0.0468 - val_loss: 0.0460
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.04601/1 [==============================] - 0s 50ms/step - loss: 0.0460 - val_loss: 0.0452
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.04521/1 [==============================] - 0s 72ms/step - loss: 0.0452 - val_loss: 0.0445
#> 1/2 [==============>...............] - ETA: 0s2/2 [==============================] - 0s 1ms/step
#> [1] "[2024-10-22 20:03:10] MXE event encoding Finish."
#> [1] "[2024-10-22 20:03:10] RI event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02571/1 [==============================] - 1s 1s/step - loss: 1.0257 - val_loss: 1.0193
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01931/1 [==============================] - 0s 60ms/step - loss: 1.0193 - val_loss: 1.0164
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01641/1 [==============================] - 0s 61ms/step - loss: 1.0164 - val_loss: 1.0147
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01471/1 [==============================] - 0s 79ms/step - loss: 1.0147 - val_loss: 1.0133
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01331/1 [==============================] - 0s 78ms/step - loss: 1.0133 - val_loss: 1.0118
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01181/1 [==============================] - 0s 76ms/step - loss: 1.0118 - val_loss: 1.0099
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00991/1 [==============================] - 0s 71ms/step - loss: 1.0099 - val_loss: 1.0071
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00711/1 [==============================] - 0s 56ms/step - loss: 1.0071 - val_loss: 1.0026
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00261/1 [==============================] - 0s 58ms/step - loss: 1.0026 - val_loss: 0.9956
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99561/1 [==============================] - 0s 59ms/step - loss: 0.9956 - val_loss: 0.9857
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98571/1 [==============================] - 0s 56ms/step - loss: 0.9857 - val_loss: 0.9735
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97351/1 [==============================] - 0s 65ms/step - loss: 0.9735 - val_loss: 0.9617
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96171/1 [==============================] - 0s 59ms/step - loss: 0.9617 - val_loss: 0.9511
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95111/1 [==============================] - 0s 83ms/step - loss: 0.9511 - val_loss: 0.9407
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94071/1 [==============================] - 0s 74ms/step - loss: 0.9407 - val_loss: 0.9287
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92871/1 [==============================] - 0s 60ms/step - loss: 0.9287 - val_loss: 0.9172
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91721/1 [==============================] - 0s 65ms/step - loss: 0.9172 - val_loss: 0.9078
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90781/1 [==============================] - 0s 69ms/step - loss: 0.9078 - val_loss: 0.9001
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90011/1 [==============================] - 0s 56ms/step - loss: 0.9001 - val_loss: 0.8935
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89351/1 [==============================] - 0s 76ms/step - loss: 0.8935 - val_loss: 0.8875
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88751/1 [==============================] - 0s 66ms/step - loss: 0.8875 - val_loss: 0.8826
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88261/1 [==============================] - 0s 56ms/step - loss: 0.8826 - val_loss: 0.8786
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87861/1 [==============================] - 0s 79ms/step - loss: 0.8786 - val_loss: 0.8744
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87441/1 [==============================] - 0s 68ms/step - loss: 0.8744 - val_loss: 0.8705
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87051/1 [==============================] - 0s 64ms/step - loss: 0.8705 - val_loss: 0.8676
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86761/1 [==============================] - 0s 75ms/step - loss: 0.8676 - val_loss: 0.8642
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86421/1 [==============================] - 0s 75ms/step - loss: 0.8642 - val_loss: 0.8608
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86081/1 [==============================] - 0s 78ms/step - loss: 0.8608 - val_loss: 0.8580
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85801/1 [==============================] - 0s 83ms/step - loss: 0.8580 - val_loss: 0.8545
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85451/1 [==============================] - 0s 75ms/step - loss: 0.8545 - val_loss: 0.8509
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85091/1 [==============================] - 0s 59ms/step - loss: 0.8509 - val_loss: 0.8477
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84771/1 [==============================] - 0s 57ms/step - loss: 0.8477 - val_loss: 0.8438
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84381/1 [==============================] - 0s 55ms/step - loss: 0.8438 - val_loss: 0.8403
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84031/1 [==============================] - 0s 64ms/step - loss: 0.8403 - val_loss: 0.8364
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83641/1 [==============================] - 0s 76ms/step - loss: 0.8364 - val_loss: 0.8325
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83251/1 [==============================] - 0s 62ms/step - loss: 0.8325 - val_loss: 0.8287
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82871/1 [==============================] - 0s 71ms/step - loss: 0.8287 - val_loss: 0.8245
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82451/1 [==============================] - 0s 60ms/step - loss: 0.8245 - val_loss: 0.8206
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82061/1 [==============================] - 0s 82ms/step - loss: 0.8206 - val_loss: 0.8164
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81641/1 [==============================] - 0s 75ms/step - loss: 0.8164 - val_loss: 0.8122
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81221/1 [==============================] - 0s 78ms/step - loss: 0.8122 - val_loss: 0.8080
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80801/1 [==============================] - 0s 81ms/step - loss: 0.8080 - val_loss: 0.8038
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80381/1 [==============================] - 0s 80ms/step - loss: 0.8038 - val_loss: 0.7995
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79951/1 [==============================] - 0s 80ms/step - loss: 0.7995 - val_loss: 0.7953
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79531/1 [==============================] - 0s 81ms/step - loss: 0.7953 - val_loss: 0.7908
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79081/1 [==============================] - 0s 80ms/step - loss: 0.7908 - val_loss: 0.7866
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78661/1 [==============================] - 0s 61ms/step - loss: 0.7866 - val_loss: 0.7823
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78231/1 [==============================] - 0s 58ms/step - loss: 0.7823 - val_loss: 0.7781
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77811/1 [==============================] - 0s 63ms/step - loss: 0.7781 - val_loss: 0.7739
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77391/1 [==============================] - 0s 77ms/step - loss: 0.7739 - val_loss: 0.7697
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76971/1 [==============================] - 0s 74ms/step - loss: 0.7697 - val_loss: 0.7656
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76561/1 [==============================] - 0s 61ms/step - loss: 0.7656 - val_loss: 0.7615
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.76151/1 [==============================] - 0s 65ms/step - loss: 0.7615 - val_loss: 0.7574
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75741/1 [==============================] - 0s 66ms/step - loss: 0.7574 - val_loss: 0.7534
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.75341/1 [==============================] - 0s 76ms/step - loss: 0.7534 - val_loss: 0.7495
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74951/1 [==============================] - 0s 63ms/step - loss: 0.7495 - val_loss: 0.7456
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74561/1 [==============================] - 0s 60ms/step - loss: 0.7456 - val_loss: 0.7418
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.74181/1 [==============================] - 0s 59ms/step - loss: 0.7418 - val_loss: 0.7380
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73801/1 [==============================] - 0s 67ms/step - loss: 0.7380 - val_loss: 0.7342
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73421/1 [==============================] - 0s 75ms/step - loss: 0.7342 - val_loss: 0.7305
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.73051/1 [==============================] - 0s 80ms/step - loss: 0.7305 - val_loss: 0.7268
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72681/1 [==============================] - 0s 78ms/step - loss: 0.7268 - val_loss: 0.7231
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.72311/1 [==============================] - 0s 84ms/step - loss: 0.7231 - val_loss: 0.7195
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71951/1 [==============================] - 0s 59ms/step - loss: 0.7195 - val_loss: 0.7160
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71601/1 [==============================] - 0s 79ms/step - loss: 0.7160 - val_loss: 0.7128
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.71281/1 [==============================] - 0s 60ms/step - loss: 0.7128 - val_loss: 0.7099
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70991/1 [==============================] - 0s 57ms/step - loss: 0.7099 - val_loss: 0.7060
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70601/1 [==============================] - 0s 57ms/step - loss: 0.7060 - val_loss: 0.7018
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.70181/1 [==============================] - 0s 58ms/step - loss: 0.7018 - val_loss: 0.6987
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69871/1 [==============================] - 0s 65ms/step - loss: 0.6987 - val_loss: 0.6955
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69551/1 [==============================] - 0s 55ms/step - loss: 0.6955 - val_loss: 0.6916
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.69161/1 [==============================] - 0s 56ms/step - loss: 0.6916 - val_loss: 0.6884
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68841/1 [==============================] - 0s 57ms/step - loss: 0.6884 - val_loss: 0.6852
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68521/1 [==============================] - 0s 57ms/step - loss: 0.6852 - val_loss: 0.6814
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.68141/1 [==============================] - 0s 79ms/step - loss: 0.6814 - val_loss: 0.6782
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67821/1 [==============================] - 0s 60ms/step - loss: 0.6782 - val_loss: 0.6749
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67491/1 [==============================] - 0s 58ms/step - loss: 0.6749 - val_loss: 0.6713
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.67131/1 [==============================] - 0s 58ms/step - loss: 0.6713 - val_loss: 0.6680
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66801/1 [==============================] - 0s 80ms/step - loss: 0.6680 - val_loss: 0.6647
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66471/1 [==============================] - 0s 73ms/step - loss: 0.6647 - val_loss: 0.6612
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.66121/1 [==============================] - 0s 57ms/step - loss: 0.6612 - val_loss: 0.6579
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65791/1 [==============================] - 0s 57ms/step - loss: 0.6579 - val_loss: 0.6547
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65471/1 [==============================] - 0s 65ms/step - loss: 0.6547 - val_loss: 0.6512
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.65121/1 [==============================] - 0s 69ms/step - loss: 0.6512 - val_loss: 0.6479
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64791/1 [==============================] - 0s 57ms/step - loss: 0.6479 - val_loss: 0.6449
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64491/1 [==============================] - 0s 76ms/step - loss: 0.6449 - val_loss: 0.6420
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.64201/1 [==============================] - 0s 64ms/step - loss: 0.6420 - val_loss: 0.6390
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63901/1 [==============================] - 0s 57ms/step - loss: 0.6390 - val_loss: 0.6352
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63521/1 [==============================] - 0s 73ms/step - loss: 0.6352 - val_loss: 0.6315
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.63151/1 [==============================] - 0s 62ms/step - loss: 0.6315 - val_loss: 0.6287
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62871/1 [==============================] - 0s 81ms/step - loss: 0.6287 - val_loss: 0.6256
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62561/1 [==============================] - 0s 75ms/step - loss: 0.6256 - val_loss: 0.6220
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.62201/1 [==============================] - 0s 72ms/step - loss: 0.6220 - val_loss: 0.6189
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61891/1 [==============================] - 0s 70ms/step - loss: 0.6189 - val_loss: 0.6160
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61601/1 [==============================] - 0s 74ms/step - loss: 0.6160 - val_loss: 0.6127
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.61271/1 [==============================] - 0s 69ms/step - loss: 0.6127 - val_loss: 0.6096
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60961/1 [==============================] - 0s 59ms/step - loss: 0.6096 - val_loss: 0.6065
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60651/1 [==============================] - 0s 75ms/step - loss: 0.6065 - val_loss: 0.6033
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60331/1 [==============================] - 0s 58ms/step - loss: 0.6033 - val_loss: 0.6002
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.60021/1 [==============================] - 0s 60ms/step - loss: 0.6002 - val_loss: 0.5973
#>  1/12 [=>............................] - ETA: 0s12/12 [==============================] - 0s 1ms/step
#> [1] "[2024-10-22 20:03:19] RI event encoding Finish."
#> [1] "[2024-10-22 20:03:19] SE event encoding..."
#> Epoch 1/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.05271/1 [==============================] - 1s 1s/step - loss: 1.0527 - val_loss: 1.0300
#> Epoch 2/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.03001/1 [==============================] - 0s 73ms/step - loss: 1.0300 - val_loss: 1.0217
#> Epoch 3/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.02171/1 [==============================] - 0s 84ms/step - loss: 1.0217 - val_loss: 1.0183
#> Epoch 4/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01831/1 [==============================] - 0s 70ms/step - loss: 1.0183 - val_loss: 1.0166
#> Epoch 5/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01661/1 [==============================] - 0s 71ms/step - loss: 1.0166 - val_loss: 1.0156
#> Epoch 6/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01561/1 [==============================] - 0s 72ms/step - loss: 1.0156 - val_loss: 1.0149
#> Epoch 7/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01491/1 [==============================] - 0s 94ms/step - loss: 1.0149 - val_loss: 1.0143
#> Epoch 8/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01431/1 [==============================] - 0s 88ms/step - loss: 1.0143 - val_loss: 1.0138
#> Epoch 9/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01381/1 [==============================] - 0s 74ms/step - loss: 1.0138 - val_loss: 1.0133
#> Epoch 10/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01331/1 [==============================] - 0s 70ms/step - loss: 1.0133 - val_loss: 1.0128
#> Epoch 11/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01281/1 [==============================] - 0s 93ms/step - loss: 1.0128 - val_loss: 1.0122
#> Epoch 12/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01221/1 [==============================] - 0s 79ms/step - loss: 1.0122 - val_loss: 1.0116
#> Epoch 13/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01161/1 [==============================] - 0s 91ms/step - loss: 1.0116 - val_loss: 1.0110
#> Epoch 14/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01101/1 [==============================] - 0s 72ms/step - loss: 1.0110 - val_loss: 1.0101
#> Epoch 15/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.01011/1 [==============================] - 0s 81ms/step - loss: 1.0101 - val_loss: 1.0091
#> Epoch 16/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00911/1 [==============================] - 0s 77ms/step - loss: 1.0091 - val_loss: 1.0077
#> Epoch 17/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00771/1 [==============================] - 0s 76ms/step - loss: 1.0077 - val_loss: 1.0059
#> Epoch 18/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00591/1 [==============================] - 0s 81ms/step - loss: 1.0059 - val_loss: 1.0034
#> Epoch 19/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00341/1 [==============================] - 0s 78ms/step - loss: 1.0034 - val_loss: 1.0001
#> Epoch 20/100
#> 1/1 [==============================] - ETA: 0s - loss: 1.00011/1 [==============================] - 0s 67ms/step - loss: 1.0001 - val_loss: 0.9961
#> Epoch 21/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99611/1 [==============================] - 0s 66ms/step - loss: 0.9961 - val_loss: 0.9918
#> Epoch 22/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.99181/1 [==============================] - 0s 68ms/step - loss: 0.9918 - val_loss: 0.9880
#> Epoch 23/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98801/1 [==============================] - 0s 80ms/step - loss: 0.9880 - val_loss: 0.9848
#> Epoch 24/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98481/1 [==============================] - 0s 68ms/step - loss: 0.9848 - val_loss: 0.9803
#> Epoch 25/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.98031/1 [==============================] - 0s 62ms/step - loss: 0.9803 - val_loss: 0.9750
#> Epoch 26/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97501/1 [==============================] - 0s 72ms/step - loss: 0.9750 - val_loss: 0.9707
#> Epoch 27/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.97071/1 [==============================] - 0s 71ms/step - loss: 0.9707 - val_loss: 0.9673
#> Epoch 28/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96731/1 [==============================] - 0s 66ms/step - loss: 0.9673 - val_loss: 0.9636
#> Epoch 29/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.96361/1 [==============================] - 0s 76ms/step - loss: 0.9636 - val_loss: 0.9594
#> Epoch 30/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95941/1 [==============================] - 0s 81ms/step - loss: 0.9594 - val_loss: 0.9553
#> Epoch 31/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95531/1 [==============================] - 0s 74ms/step - loss: 0.9553 - val_loss: 0.9522
#> Epoch 32/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.95221/1 [==============================] - 0s 77ms/step - loss: 0.9522 - val_loss: 0.9495
#> Epoch 33/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94951/1 [==============================] - 0s 76ms/step - loss: 0.9495 - val_loss: 0.9464
#> Epoch 34/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94641/1 [==============================] - 0s 78ms/step - loss: 0.9464 - val_loss: 0.9433
#> Epoch 35/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94331/1 [==============================] - 0s 70ms/step - loss: 0.9433 - val_loss: 0.9407
#> Epoch 36/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.94071/1 [==============================] - 0s 85ms/step - loss: 0.9407 - val_loss: 0.9383
#> Epoch 37/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93831/1 [==============================] - 0s 85ms/step - loss: 0.9383 - val_loss: 0.9355
#> Epoch 38/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93551/1 [==============================] - 0s 69ms/step - loss: 0.9355 - val_loss: 0.9326
#> Epoch 39/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93261/1 [==============================] - 0s 82ms/step - loss: 0.9326 - val_loss: 0.9301
#> Epoch 40/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.93011/1 [==============================] - 0s 86ms/step - loss: 0.9301 - val_loss: 0.9273
#> Epoch 41/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92731/1 [==============================] - 0s 69ms/step - loss: 0.9273 - val_loss: 0.9244
#> Epoch 42/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92441/1 [==============================] - 0s 81ms/step - loss: 0.9244 - val_loss: 0.9217
#> Epoch 43/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.92171/1 [==============================] - 0s 78ms/step - loss: 0.9217 - val_loss: 0.9188
#> Epoch 44/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91881/1 [==============================] - 0s 65ms/step - loss: 0.9188 - val_loss: 0.9158
#> Epoch 45/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91581/1 [==============================] - 0s 92ms/step - loss: 0.9158 - val_loss: 0.9131
#> Epoch 46/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91311/1 [==============================] - 0s 88ms/step - loss: 0.9131 - val_loss: 0.9102
#> Epoch 47/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.91021/1 [==============================] - 0s 64ms/step - loss: 0.9102 - val_loss: 0.9073
#> Epoch 48/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90731/1 [==============================] - 0s 65ms/step - loss: 0.9073 - val_loss: 0.9045
#> Epoch 49/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90451/1 [==============================] - 0s 64ms/step - loss: 0.9045 - val_loss: 0.9016
#> Epoch 50/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.90161/1 [==============================] - 0s 66ms/step - loss: 0.9016 - val_loss: 0.8988
#> Epoch 51/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89881/1 [==============================] - 0s 73ms/step - loss: 0.8988 - val_loss: 0.8960
#> Epoch 52/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89601/1 [==============================] - 0s 64ms/step - loss: 0.8960 - val_loss: 0.8931
#> Epoch 53/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89311/1 [==============================] - 0s 67ms/step - loss: 0.8931 - val_loss: 0.8903
#> Epoch 54/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.89031/1 [==============================] - 0s 66ms/step - loss: 0.8903 - val_loss: 0.8875
#> Epoch 55/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88751/1 [==============================] - 0s 82ms/step - loss: 0.8875 - val_loss: 0.8847
#> Epoch 56/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88471/1 [==============================] - 0s 73ms/step - loss: 0.8847 - val_loss: 0.8819
#> Epoch 57/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.88191/1 [==============================] - 0s 90ms/step - loss: 0.8819 - val_loss: 0.8792
#> Epoch 58/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87921/1 [==============================] - 0s 75ms/step - loss: 0.8792 - val_loss: 0.8764
#> Epoch 59/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87641/1 [==============================] - 0s 106ms/step - loss: 0.8764 - val_loss: 0.8737
#> Epoch 60/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87371/1 [==============================] - 0s 89ms/step - loss: 0.8737 - val_loss: 0.8709
#> Epoch 61/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.87091/1 [==============================] - 0s 78ms/step - loss: 0.8709 - val_loss: 0.8681
#> Epoch 62/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86811/1 [==============================] - 0s 94ms/step - loss: 0.8681 - val_loss: 0.8653
#> Epoch 63/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86531/1 [==============================] - 0s 104ms/step - loss: 0.8653 - val_loss: 0.8625
#> Epoch 64/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.86251/1 [==============================] - 0s 86ms/step - loss: 0.8625 - val_loss: 0.8597
#> Epoch 65/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85971/1 [==============================] - 0s 73ms/step - loss: 0.8597 - val_loss: 0.8569
#> Epoch 66/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85691/1 [==============================] - 0s 105ms/step - loss: 0.8569 - val_loss: 0.8541
#> Epoch 67/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85411/1 [==============================] - 0s 72ms/step - loss: 0.8541 - val_loss: 0.8512
#> Epoch 68/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.85121/1 [==============================] - 0s 66ms/step - loss: 0.8512 - val_loss: 0.8484
#> Epoch 69/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84841/1 [==============================] - 0s 68ms/step - loss: 0.8484 - val_loss: 0.8456
#> Epoch 70/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84561/1 [==============================] - 0s 99ms/step - loss: 0.8456 - val_loss: 0.8428
#> Epoch 71/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84281/1 [==============================] - 0s 88ms/step - loss: 0.8428 - val_loss: 0.8400
#> Epoch 72/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.84001/1 [==============================] - 0s 104ms/step - loss: 0.8400 - val_loss: 0.8373
#> Epoch 73/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83731/1 [==============================] - 0s 89ms/step - loss: 0.8373 - val_loss: 0.8345
#> Epoch 74/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83451/1 [==============================] - 0s 116ms/step - loss: 0.8345 - val_loss: 0.8318
#> Epoch 75/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.83181/1 [==============================] - 0s 79ms/step - loss: 0.8318 - val_loss: 0.8291
#> Epoch 76/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82911/1 [==============================] - 0s 63ms/step - loss: 0.8291 - val_loss: 0.8264
#> Epoch 77/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82641/1 [==============================] - 0s 77ms/step - loss: 0.8264 - val_loss: 0.8238
#> Epoch 78/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82381/1 [==============================] - 0s 83ms/step - loss: 0.8238 - val_loss: 0.8212
#> Epoch 79/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.82121/1 [==============================] - 0s 68ms/step - loss: 0.8212 - val_loss: 0.8186
#> Epoch 80/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81861/1 [==============================] - 0s 79ms/step - loss: 0.8186 - val_loss: 0.8160
#> Epoch 81/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81601/1 [==============================] - 0s 83ms/step - loss: 0.8160 - val_loss: 0.8135
#> Epoch 82/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81351/1 [==============================] - 0s 102ms/step - loss: 0.8135 - val_loss: 0.8110
#> Epoch 83/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.81101/1 [==============================] - 0s 90ms/step - loss: 0.8110 - val_loss: 0.8085
#> Epoch 84/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80851/1 [==============================] - 0s 105ms/step - loss: 0.8085 - val_loss: 0.8061
#> Epoch 85/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80611/1 [==============================] - 0s 86ms/step - loss: 0.8061 - val_loss: 0.8038
#> Epoch 86/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80381/1 [==============================] - 0s 93ms/step - loss: 0.8038 - val_loss: 0.8018
#> Epoch 87/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.80181/1 [==============================] - 0s 88ms/step - loss: 0.8018 - val_loss: 0.7998
#> Epoch 88/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79981/1 [==============================] - 0s 82ms/step - loss: 0.7998 - val_loss: 0.7973
#> Epoch 89/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79731/1 [==============================] - 0s 77ms/step - loss: 0.7973 - val_loss: 0.7944
#> Epoch 90/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79441/1 [==============================] - 0s 60ms/step - loss: 0.7944 - val_loss: 0.7926
#> Epoch 91/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79261/1 [==============================] - 0s 82ms/step - loss: 0.7926 - val_loss: 0.7907
#> Epoch 92/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.79071/1 [==============================] - 0s 63ms/step - loss: 0.7907 - val_loss: 0.7881
#> Epoch 93/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78811/1 [==============================] - 0s 64ms/step - loss: 0.7881 - val_loss: 0.7860
#> Epoch 94/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78601/1 [==============================] - 0s 70ms/step - loss: 0.7860 - val_loss: 0.7842
#> Epoch 95/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78421/1 [==============================] - 0s 67ms/step - loss: 0.7842 - val_loss: 0.7819
#> Epoch 96/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.78191/1 [==============================] - 0s 71ms/step - loss: 0.7819 - val_loss: 0.7796
#> Epoch 97/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77961/1 [==============================] - 0s 73ms/step - loss: 0.7796 - val_loss: 0.7779
#> Epoch 98/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77791/1 [==============================] - 0s 65ms/step - loss: 0.7779 - val_loss: 0.7758
#> Epoch 99/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77581/1 [==============================] - 0s 64ms/step - loss: 0.7758 - val_loss: 0.7736
#> Epoch 100/100
#> 1/1 [==============================] - ETA: 0s - loss: 0.77361/1 [==============================] - 0s 80ms/step - loss: 0.7736 - val_loss: 0.7718
#>  1/61 [..............................] - ETA: 4s45/61 [=====================>........] - ETA: 0s61/61 [==============================] - 0s 1ms/step
#> [1] "[2024-10-22 20:03:30] SE event encoding Finish."
#> [1] "[2024-10-22 20:03:30] step5 Calculate splicing regulation distance and Combine distance ======="
#> [1] "384 rbps are used to calculate splicing regulation information"
#> [1] "Save data"
#> [1] "Save data Finished"
#> [1] "[2024-10-22 20:04:48] step6 Calculate combined event similarity ======="
#> [1] "[2024-10-22 20:04:49] Calculate  A3SS event Similarity"
#> 
#> Attaching package: 'Matrix'
#> The following object is masked from 'package:S4Vectors':
#> 
#>     expand
#> [1] "[2024-10-22 20:04:59] Calculate A3SS event Similarity Finished"
#> [1] "[2024-10-22 20:04:59] Calculate  A5SS event Similarity"
#> [1] "[2024-10-22 20:05:08] Calculate A5SS event Similarity Finished"
#> [1] "[2024-10-22 20:05:09] Calculate  MXE event Similarity"
#> [1] "[2024-10-22 20:05:18] Calculate MXE event Similarity Finished"
#> [1] "[2024-10-22 20:05:19] Calculate  RI event Similarity"
#> [1] "[2024-10-22 20:05:28] Calculate RI event Similarity Finished"
#> [1] "[2024-10-22 20:05:29] Calculate  SE event Similarity"
#> [1] "[2024-10-22 20:05:38] Calculate SE event Similarity Finished"
#> [1] "[2024-10-22 20:05:38] Calculate event similarity Finished."
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
#> [1] "[2024-10-22 20:05:39] Get imputed result using cell similarity and event similarity."
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
#> [1] "[2024-10-22 20:05:40] Running Event_type=A3SS;cell_similarity_feature=PSI"
#> [1] "[2024-10-22 20:05:40] Save data"
#> [1] "[2024-10-22 20:05:40] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500020898.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500020898.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:05:54] Running Event_type=A3SS;cell_similarity_feature=RC"
#> [1] "[2024-10-22 20:05:54] Save data"
#> [1] "[2024-10-22 20:05:54] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500018759.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500018759.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:06:07] Running Event_type=A3SS;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-22 20:06:07] Save data"
#> [1] "[2024-10-22 20:06:07] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500008554.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500008554.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:06:21] Running Event_type=A5SS;cell_similarity_feature=PSI"
#> [1] "[2024-10-22 20:06:21] Save data"
#> [1] "[2024-10-22 20:06:21] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500001430.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500001430.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:06:34] Running Event_type=A5SS;cell_similarity_feature=RC"
#> [1] "[2024-10-22 20:06:34] Save data"
#> [1] "[2024-10-22 20:06:34] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500009030.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500009030.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:06:47] Running Event_type=A5SS;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-22 20:06:47] Save data"
#> [1] "[2024-10-22 20:06:47] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500002168.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500002168.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:07:01] Running Event_type=MXE;cell_similarity_feature=PSI"
#> [1] "[2024-10-22 20:07:01] Save data"
#> [1] "[2024-10-22 20:07:01] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500003874.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500003874.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:07:14] Running Event_type=MXE;cell_similarity_feature=RC"
#> [1] "[2024-10-22 20:07:14] Save data"
#> [1] "[2024-10-22 20:07:14] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-499990105.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-499990105.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:07:27] Running Event_type=MXE;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-22 20:07:27] Save data"
#> [1] "[2024-10-22 20:07:27] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500020646.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500020646.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:07:40] Running Event_type=RI;cell_similarity_feature=PSI"
#> [1] "[2024-10-22 20:07:40] Save data"
#> [1] "[2024-10-22 20:07:41] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-499996850.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-499996850.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:07:54] Running Event_type=RI;cell_similarity_feature=RC"
#> [1] "[2024-10-22 20:07:54] Save data"
#> [1] "[2024-10-22 20:07:54] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500019703.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500019703.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:08:07] Running Event_type=RI;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-22 20:08:07] Save data"
#> [1] "[2024-10-22 20:08:07] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-499992560.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-499992560.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:08:20] Running Event_type=SE;cell_similarity_feature=PSI"
#> [1] "[2024-10-22 20:08:20] Save data"
#> [1] "[2024-10-22 20:08:20] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500005041.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500005041.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:08:33] Running Event_type=SE;cell_similarity_feature=RC"
#> [1] "[2024-10-22 20:08:33] Save data"
#> [1] "[2024-10-22 20:08:33] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-499978980.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-499978980.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:08:46] Running Event_type=SE;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-10-22 20:08:46] Save data"
#> [1] "[2024-10-22 20:08:46] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_3871154-500012937.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_3871154-500012937.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-10-22 20:09:00] Get imputed result using cell similarity and event similarity Finish."
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
#>   ..$ PSI_PSI    : num [1:2550, 1:15] 0.696 0.634 0.586 0.501 0.581 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ RC_PSI     : num [1:2550, 1:15] 0.658 0.611 0.575 0.496 0.565 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ EXP_RBP_PSI: num [1:2550, 1:15] 0.616 0.603 0.573 0.527 0.584 ...
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
#> [1] "[2024-10-22 20:09:00] Combine imputed psi."
#> [1] "[2024-10-22 20:09:00] Loading data..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_seperated_500012652.rds"
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
#> [1] "[2024-10-22 20:09:04] Combine imputed psi Finish."
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
#> [1] "[2024-10-22 20:09:04] Counting reads of A3SS events..."
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
#> [1] "[2024-10-22 20:09:06] Counting reads of A3SS events Finish."
#> [1] "[2024-10-22 20:09:06] Counting reads of A5SS events..."
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
#> [1] "[2024-10-22 20:09:09] Counting reads of A5SS events Finish."
#> [1] "[2024-10-22 20:09:09] Counting reads of MXE events..."
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
#> [1] "[2024-10-22 20:09:10] Counting reads of MXE events Finish."
#> [1] "[2024-10-22 20:09:10] Counting reads of SE events..."
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
#> [1] "[2024-10-22 20:09:14] Counting reads of SE events Finish."
ftpsi.path = getFtRawPSI(paras)
#> [1] "Checking raw reads..."
#> [1] "Loading splicing events for classifer fine tune..."
#> [1] "[2024-10-22 20:09:14] Calculating PSI value of A3SS events..."
#> [1] "[2024-10-22 20:09:14] Calculating PSI value of A3SS events Finish."
#> [1] "[2024-10-22 20:09:14] Calculating PSI value of A5SS events..."
#> [1] "[2024-10-22 20:09:14] Calculating PSI value of A5SS events Finish."
#> [1] "[2024-10-22 20:09:14] Calculating PSI value of MXE events..."
#> [1] "[2024-10-22 20:09:14] Calculating PSI value of MXE events Finish."
#> [1] "[2024-10-22 20:09:14] Calculating PSI value of SE events..."
#> [1] "[2024-10-22 20:09:14] Calculating PSI value of SE events Finish."
ftrds.path = mergeFtSplicingValue(paras)
ftmodel.path = FtClassifier(paras)
#> [1] "Reading true Ft PSI..."
#> [1] "Loading Pre-training classifer..."
#> [1] "[2024-10-22 20:09:14] Classifer fine tune"
#> [1] "[2024-10-22 20:09:14] Processing raw Ft data..."
#> [1] "Checking data..."
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=PSI;RC;EXP_RBP  checked"
#> [1] "Calculating Classifier Features..."
#> [1] "[2024-10-22 20:09:15] Save data"
#> [1] "[2024-10-22 20:09:15] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-499998203.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-499998203.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:09:38] Save data"
#> [1] "[2024-10-22 20:09:38] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-500013636.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-500013636.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:10:01] Save data"
#> [1] "[2024-10-22 20:10:01] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-500011055.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-500011055.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:10:25] Save data"
#> [1] "[2024-10-22 20:10:25] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-499990816.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-499990816.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:10:48] Save data"
#> [1] "[2024-10-22 20:10:48] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-500011480.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-500011480.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:11:11] Save data"
#> [1] "[2024-10-22 20:11:11] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-500010175.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-500010175.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:11:35] Save data"
#> [1] "[2024-10-22 20:11:35] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-499980377.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-499980377.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:11:59] Save data"
#> [1] "[2024-10-22 20:11:59] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-500002586.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-500002586.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:12:24] Save data"
#> [1] "[2024-10-22 20:12:24] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-499990265.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-499990265.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:12:48] Save data"
#> [1] "[2024-10-22 20:12:48] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-499994241.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-499994241.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:13:12] Save data"
#> [1] "[2024-10-22 20:13:12] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-500010052.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-500010052.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:13:36] Save data"
#> [1] "[2024-10-22 20:13:36] Save data Finished"
#> [1] "bash /tmp/Rtmp3pLKHL/temp_libpath3afa5936e27ff7/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_3871154-499984500.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_3871154-499984500.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-10-22 20:13:59]  Model training;similarity_type=PSI"
#> [1] "[2024-10-22 20:14:00]  Model training;similarity_type=RC"
#> [1] "[2024-10-22 20:14:01]  Model training;similarity_type=EXP_RBP"
#> [1] "[2024-10-22 20:14:02] Classifer fine tune Finish."
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
ImputedFt.data.final.path = Estimation(paras,rds_imputed_file = Imputed.data.path)
#> [1] "[2024-10-22 20:14:02] Combine imputed psi."
#> [1] "[2024-10-22 20:14:02] Loading data..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_seperated_500012652.rds"
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
#> [1] "[2024-10-22 20:14:06] Combine imputed psi Finish."

print(ImputedFt.data.final.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_combined_499969497.rds"
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
