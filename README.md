![alt text](https://github.com/lvxuan12/SCSES/blob/main/png/SCSES_logo.png)
Menu
================
- [SCSES](#scses)
  - [Hardware requirements](#hardware-requirements)
  - [Software requirements](#software-requirements)
    - [OS Requirements](#os-requirements)
  - [Installation](#installation)
    - [Installation with docker file](#Installation-with-docker-file)
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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCSES

<!-- badges: start -->
<!-- badges: end -->

Single-cell Splicing Estimation based on Network Diffusion

## Hardware requirements

`SCSES` package requires only a standard computer with enough RAM to
support the in-memory operations.

## Software requirements

### OS Requirements

This package is supported for Linux. The package has been tested on
Debian-11.21.

## Installation
### Installation with docker file
SCSES provides a Docker-based installation method to simplify the setup of all dependencies and requirements. Please follow the steps below to build the Docker image and start a container to use SCSES:
#### 1. Install Docker Client: 
Please install the [Docker client](https://www.docker.com/products/docker-desktop) on the host machine.
#### 2. Download Dockerfile
The Dockerfile of SCSES can be downloaded from: https://github.com/lvxuan12/SCSES/blob/main/SCSES.dockerfile.
#### 3.	Build Docker Image
You can build the SCSES Docker image using the command:
``` bash
docker build -t scses -f SCSES.dockerfile .
```

#### 4. Create Docker Container
After building the image, create a Docker container with the following command:
``` bash
docker run -d -p [exported port]:8787 -e PASSWORD=[user password] -v [local directory]:/data --name test scses
```
``[exported port]``: The port on the host machine to access the container.  
``[user password]``: A user-defined password for logging into the RStudio server.  
``[local directory]``: A local directory mapped to the container for data storage and sharing.  
#### 5.	Access RStudio Server
Now, you can access the RStudio server by opening a web browser and navigating to ``[host IP]:[exported port]``. Use the default username ``rstudio`` and the user-defined password to log in.  
In this pre-configured RStudio server environment, SCSES and all its dependencies are correctly installed and ready for use.  
Please refer to [Getting started][#getting-started] to start the first experience with SCSES. Enjoy!
### Installation of dependencies and requirements

We recommend a new **conda** environment to install SCSES:

``` bash
conda create -n SCSES_test python=3.11
conda activate SCSES_test
```

To use SCSES, you will need to install R, Python, Matlab Compiler
Runtime(v9.13), and Java(v17.0.10).

The MCR is quite large, so downloading may take some time.

``` bash
## install R in conda environment
conda install -c conda-forge r-base=4.3.1
## install MCR
mkdir /path/to/MCR && \
cd /path/to/MCR && \
wget https://ssd.mathworks.com/supportfiles/downloads/R2022b/Release/10/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2022b_Update_10_glnxa64.zip && \
unzip -q MATLAB_Runtime_R2022b_Update_10_glnxa64.zip && \
./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent
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
cd rmats_turbo_v4_3_0/
pip install Cython
./build_rmats
export PATH=/path/to/rmats_turbo_v4_3_0/:$PATH
```

##### 2.2 [MAJIQ](https://biociphers.bitbucket.io/majiq-docs/index.html)

``` bash
conda create -n MAJIQ python=3.11
conda activate MAJIQ

# install HTSlib from <https://www.htslib.org/download/>
# change to where library/include directories are installed

export HTSLIB_LIBRARY_DIR=/path/to/htslib-1.15.1/lib/
export HTSLIB_INCLUDE_DIR=/path/to/htslib-1.15.1/include/

pip install git+https://bitbucket.org/biociphers/majiq_academic.git
export MAJIQ_LICENSE_FILE=/path/to/majiq_license_academic_official.lic
```

**NOTE:** MAJIQ will not function without providing the license file.

##### 2.3 [IRFinder](https://github.com/dgaolab/IRFinder)

``` bash
wget https://github.com/RitchieLabIGH/IRFinder/archive/refs/tags/v2.0.1.tar.gz
tar -zxvf v2.0.1.tar.gz
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
install.packages("remotes")
remotes::install_version("Matrix", version = "1.6-5")
install.packages("curl",configure.vars='LIB_DIR=/usr/lib/x86_64-linux-gnu/pkgconfig/')
options(download.file.method = "wget", times=100)
remotes::install_github("lvxuan12/SCSES")
```

#### Tips for some Installation errors

##### 1. cannot find fftw.h

``` bash
conda install conda-forge::fftw
# Configure FFTW
export FFTW_CFLAGS=" -I/path/to/miniconda3/envs/scses/include/"
export FFTW_LIBS=" -L/path/to/miniconda3/envs/scses/lib -lfftw3"
```

##### 2. cannot find -lxml2

``` bash
conda install conda-forge::libxml2
```

##### 3. cannot find -lsz

``` bash
ln -s /usr/lib/x86_64-linux-gnu/libsz.so /path/to/miniconda/envs/SCSES_test/lib/libsz.so
```

##### 4. error: ‘::timespec_get’ has not been declared

``` bash
conda upgrade -c conda-forge --all 
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
#> Bioconductor version '3.16' is out-of-date; the current release version '3.20'
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

Move test bam files to directory `/path/to/bams/`.

Move other input data to directory `/path/to/refgenome/`.

``` bash
# Example:
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
## Loading packages
library(SCSES)
#paras_file: path to configure file generated in the previous step
paras = readSCSESconfig(paras_file)
```

The `cell_line.json` file is an example configuration file for test data
which can be downloaded previously or load from SCSES package
[here](https://github.com/lvxuan12/SCSES/blob/main/analysis/cell_line.json).

For real dataset, users can modify this file to fit their input and
software environment or use `createConfigshiny` function to create a new
configuration file.

``` r
## Loading packages
library(SCSES)
#paras_file: path to configure file generated in the previous step
print(system.file("analysis", package = "SCSES"))
#> [1] "/tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/analysis"
paras = readSCSESconfig(paste0(system.file("analysis", package = "SCSES"),'/cell_line.json'))
#paras = readSCSESconfig("/disk/lvxuan/Single-Splicing/src/package/SCSES/analysis/cell_line.json")
names(paras)
#> [1] "DataSet" "Basic"   "Task"
print(paras$DataSet)
#> [1] "cell_line"
print(paras$Basic$bam_path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam"
# all outputs will be saved in work_path
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
rds.path = getEXPmatrix(paras)
```

``` r
featurecounts.path = getGeneExpression(paras) 
#> [1] "[2024-11-17 15:10:36] Detect gene expression: bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/shell/run_featurecounts.sh /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/expr/ /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/test.fa /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/test.gtf /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam 20 cell_line paired /disk/software/subread-2.0.6-source/bin/featureCounts >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/runfeatureCounts.log 2>&1"
#> [1] "[2024-11-17 15:10:54] Detect gene expression Finish."
rds.path = getEXPmatrix(paras)
print(rds.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds/"
list.files(rds.path)
#> [1] "count_norm.rds" "count.rds"      "event.rds"      "psi.rds"       
#> [5] "rc.rds"         "TPM.rds"
tpm = readRDS(paste0(rds.path,'/count_norm.rds'))
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

After this step, directories `expr` and `rds` will be created.

#### Normalized UMI count matrix (for UMI dataset)

You can use `get10XEXPmatrix` to generate Normalized UMI count matrix
from 10X CellRanger hdf5 file, which will save normalized UMI count to
`work_path/rds/`.

``` r
# install.packages('Seurat')
# In expr_path, there are different subdirectories, each named after a sample name.
rds.path = get10XEXPmatrix(paras,expr_path,sample_name)
```

### Step3. Detect splicing events

To define a global set of all splicing events, SCSES firstly merges all
bam files from every single cell to construct a pseudo-bulk bam file,
and identifies all types of splicing events by conventional algorithms.

#### for smart-seq2 dataset

``` r
pseudobulk.path = createPseudobulk(paras)
event.path = detectEvents(paras,star_ref_path)
```

``` r
pseudobulk.path = createPseudobulk(paras)
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data/"
#> [1] "[2024-11-17 15:10:55] Creating Pseudobulk directory..."
#> Warning in dir.create(path = pseudobulk.path, recursive = T):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data' already exists
#> [1] "Pseudobulk bam file exists."
#> [1] "[2024-11-17 15:10:55] Merge Bam Files: /disk/software/samtools/bin/samtools merge -f -@ 20 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data//all.bam /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/bam/*.bam --no-PG"
#> [1] "[2024-11-17 15:11:37] Merge Bam Files Finish."
#> [1] "[2024-11-17 15:11:37] Bam File Index: /disk/software/samtools/bin/samtools index -@ 20 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data//all.bam"
#> [1] "[2024-11-17 15:11:43] Bam File Index Finish."
print(pseudobulk.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/data/"
list.files(pseudobulk.path)
#> [1] "all.bam"     "all.bam.bai"

#if you meet:
#irfinder: error while loading shared libraries: libboost_iostreams.so.1.71.0: cannot open shared object file: No such file or directory
#libboost_iostreams.so.1.71.0 exists in /disk/lvxuan/lib
old.ld=Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH = paste0("/disk/lvxuan/lib:", old.ld))
# you can Provide STAR reference to speed up the function
event.path = detectEvents(paras,star_ref_path="/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/refgenome/STAR_Reference")
#> [1] "Checking cells..."
#> [1] "Checking events..."
#> [1] "event_type=SE;RI;A3SS;A5SS;MXE  checked"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/events/"
#> [1] "[2024-11-17 15:11:43] Creating events directory..."
#> Warning in dir.create(path = work_path):
#> '/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/events' already
#> exists
#> [1] "[2024-11-17 15:11:43] Generating SE event id"
#> arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
#> [1] "[2024-11-17 15:11:43] Generating SE event id Finish."
#> [1] "[2024-11-17 15:11:43] Generating RI event id"
#> arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
#> [1] "[2024-11-17 15:11:59] Generating RI event id Finish."
#> [1] "[2024-11-17 15:11:59] Generating A3SS event id"
#> [1] "[2024-11-17 15:11:59] Generating A3SS event id Finish."
#> [1] "[2024-11-17 15:11:59] Generating A5SS event id"
#> [1] "[2024-11-17 15:11:59] Generating A5SS event id Finish."
#> [1] "[2024-11-17 15:11:59] Generating MXE event id"
#> [1] "[2024-11-17 15:11:59] Generating MXE event id Finish."
#> [1] "Total splicing event: SE=5764 RI=582 A3SS=237 A5SS=231 MXE=349"
print(event.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/events/"
list.files(event.path)
#>  [1] "A3SS.txt"            "A5SS.txt"            "IRFinder"           
#>  [4] "java_getRC_A3SS.log" "java_getRC_A5SS.log" "java_getRC_MXE.log" 
#>  [7] "java_getRC_RI.log"   "java_getRC_SE.log"   "majiq"              
#> [10] "MXE.txt"             "RI.txt"              "rMats"              
#> [13] "SE.txt"
## SE events
se.event=readLines(paste0(event.path,'SE.txt'))
print(nrow(se.event))
#> NULL
print(se.event[1])
#> [1] "isoform1=exon:chr5:180219870-180220097:-@junction:chr5:180220098-180229679:-@exon:chr5:180229680-180229803:-|isoform2=junction:chr5:180220098-180222655:-@exon:chr5:180222656-180222875:-@junction:chr5:180222876-180229679:-|MGAT1|SE"
```

Different types of splicing events will be saved to `work_path/events/`,
separately.

#### for UMI dataset

SCSES requires single cell bam files being saved in a directory. For
UMI-based dataset using CellRanger for data process, the function
`split10XBAM` can be used to get single cell bam files for each sample.

``` r
# CellRanger_path: directory to CellRanger output
# out_path: directory to save single cell bam
# core: the number of threads

splitbam.path = split10XBAM(CellRanger_path,out_path,core)

# path to single cell bam files should be added to bam_path in configure file
paras$Basic$bam_path = splitbam.path
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
rawpsi.path = getRawPSI(paras)
rawrds.path = mergeSplicingValue(paras)
processed.data.path = preprocessEvent(paras)
```

``` r
rawrc.path = getRawRC(paras)
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "Checking cells..."
#> [1] "15 cells are considered."
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value/"
#> [1] "Splicing event types: A3SS A5SS MXE RI SE"
#> [1] "[2024-11-17 15:11:59] Counting reads of A3SS events..."
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
#> [1] "[2024-11-17 15:12:02] Counting reads of A3SS events Finish."
#> [1] "[2024-11-17 15:12:02] Counting reads of A5SS events..."
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
#> [1] "[2024-11-17 15:12:05] Counting reads of A5SS events Finish."
#> [1] "[2024-11-17 15:12:05] Counting reads of MXE events..."
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
#> [1] "[2024-11-17 15:12:10] Counting reads of MXE events Finish."
#> [1] "[2024-11-17 15:12:10] Counting reads of RI events..."
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
#> [1] "[2024-11-17 15:12:19] Counting reads of RI events Finish."
#> [1] "[2024-11-17 15:12:19] Counting reads of SE events..."
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
#> [1] "[2024-11-17 15:12:48] Counting reads of SE events Finish."
rawpsi.path = getRawPSI(paras)
#> [1] "Checking raw reads..."
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of A3SS events..."
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of A3SS events Finish."
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of A5SS events..."
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of A5SS events Finish."
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of MXE events..."
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of MXE events Finish."
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of RI events..."
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of RI events Finish."
#> [1] "[2024-11-17 15:12:48] Calculating PSI value of SE events..."
#> [1] "[2024-11-17 15:12:49] Calculating PSI value of SE events Finish."
rawrds.path = mergeSplicingValue(paras)
processed.data.path = preprocessEvent(paras)
#> [1] "[2024-11-17 15:12:49] Processing raw data..."
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
#> [1] "[2024-11-17 15:12:50] Successfully processed data."
print(rawrds.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds/"
print(processed.data.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds_processed/"

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
## merged raw data
psi = readRDS(paste0(rawrds.path,'/psi.rds'))
print(psi[1:3,1:3])
#>                                                                                                                                                                                                                                        SRR11826368.bam
#> isoform1=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178086787:+@exon:chr2:178086788-178088686:+|isoform2=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178084135:+@exon:chr2:178084136-178088686:+|HNRNPA3|A3SS       1.0000000
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS         0.8888889
#> isoform1=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041195:+@exon:chr3:184041196-184041381:+|isoform2=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041192:+@exon:chr3:184041193-184041381:+|EIF4G1|A3SS        0.4615385
#>                                                                                                                                                                                                                                        SRR11826371.bam
#> isoform1=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178086787:+@exon:chr2:178086788-178088686:+|isoform2=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178084135:+@exon:chr2:178084136-178088686:+|HNRNPA3|A3SS       1.0000000
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS         1.0000000
#> isoform1=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041195:+@exon:chr3:184041196-184041381:+|isoform2=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041192:+@exon:chr3:184041193-184041381:+|EIF4G1|A3SS        0.2727273
#>                                                                                                                                                                                                                                        SRR11826409.bam
#> isoform1=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178086787:+@exon:chr2:178086788-178088686:+|isoform2=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178084135:+@exon:chr2:178084136-178088686:+|HNRNPA3|A3SS       1.0000000
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS         0.8888889
#> isoform1=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041195:+@exon:chr3:184041196-184041381:+|isoform2=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041192:+@exon:chr3:184041193-184041381:+|EIF4G1|A3SS        0.2500000
print(dim(psi))
#> [1] 7149   16

rc = readRDS(paste0(rawrds.path,'/rc.rds'))
print(rc[1:3,1:3])
#>                                     SRR11826368.bam SRR11826371.bam
#> junction:chr2:178084043-178086787:+               0               0
#> junction:chr4:146019572-146025540:+               2               0
#> junction:chr3:184041030-184041195:+              14               8
#>                                     SRR11826409.bam
#> junction:chr2:178084043-178086787:+               0
#> junction:chr4:146019572-146025540:+               1
#> junction:chr3:184041030-184041195:+              21
print(dim(rc))
#> [1] 16049    16

event = readRDS(paste0(rawrds.path,'/event.rds'))
print(event[1:3,])
#>                                                                                                                                                                                                                                    event
#> 1 isoform1=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178086787:+@exon:chr2:178086788-178088686:+|isoform2=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178084135:+@exon:chr2:178084136-178088686:+|HNRNPA3|A3SS
#> 2   isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS
#> 3  isoform1=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041195:+@exon:chr3:184041196-184041381:+|isoform2=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041192:+@exon:chr3:184041193-184041381:+|EIF4G1|A3SS
#>                            exclusion1 exclusion2
#> 1 junction:chr2:178084043-178086787:+       <NA>
#> 2 junction:chr4:146019572-146025540:+       <NA>
#> 3 junction:chr3:184041030-184041195:+       <NA>
#>                            retention1 retention2 type
#> 1 junction:chr2:178084043-178084135:+       <NA> A3SS
#> 2 junction:chr4:146019572-146025537:+       <NA> A3SS
#> 3 junction:chr3:184041030-184041192:+       <NA> A3SS

## processed data
psi_processed = readRDS(paste0(processed.data.path,'/psi.rds'))
rc_processed = readRDS(paste0(processed.data.path,'/rc.rds'))
print(dim(psi_processed))
#> [1] 2550   15
print(dim(rc_processed))
#> [1] 6087   15
```

Raw read count matrix, PSI matrix, and event annotation will be saved to
`work_path/rds/`.Then, data after quality control process will be saved
to `work_path/rds_processed/`, which will be used for **subsequent
calculations**.

### Step5. Constructs similarity networks

To overcome the high dropout rate and limited read coverage of scRNA-seq
techniques, SCSES constructs cell similarity and event similarity
networks by K-nearest neighbor algorithm (KNN) to learn information from
similar cells/events.

``` r
cellnet.path = getCellSimilarity(paras)
eventnet.path = getEventSimilarity(paras)
```

``` r
cellnet.path = getCellSimilarity(paras)
#> [1] "[2024-11-17 15:12:51] Calculate cell similarity..."
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
#> [1] "[2024-11-17 15:12:51] Computing cell similarity based on PSI"
#> [1] "Calculate similarity among 15 cells."
#> Delta: [0.17081888 0.02261298]
#> [1] "[2024-11-17 15:12:59] Computing cell similarity based on RC"
#> [1] "Calculate similarity among 15 cells."
#> Delta: [0.13512509 0.02119545]
#> [1] "[2024-11-17 15:12:59] Computing cell similarity based on EXP_RBP"
#> [1] "Calculate similarity among 15 cells."
#> [1] "The number of features is greater than the number of rows in the input data."
#> [1] "Total 384 features will be used"
#> Delta: [0.17068148 0.03534559]
#> [1] "[2024-11-17 15:12:59] Calculate cell similarity Finish."
eventnet.path = getEventSimilarity(paras)
#> [1] "[2024-11-17 15:12:59] Calculate events KNN..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/rds_processed/"
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation/event_similarity/"
#> [1] "alpha_event=0.8  checked"
#> [1] "kevent=5  checked"
#> [1] "decay_event=0.05  checked"
#> [1] "Checking data..."
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "[2024-11-17 15:13:01] Calculate events feature..."
#> [1] "[2024-11-17 15:13:01] step1 Creating BSgenome for hg19======="
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
#> [1] "[2024-11-17 15:13:02] step2 Extracting features ======="
#> [1] "[2024-11-17 15:13:02] Extracting A3SS features..."
#> [1] "[2024-11-17 15:13:02] Loading events..."
#> [1] "[2024-11-17 15:13:12] Parsing events region..."
#> [1] "[2024-11-17 15:13:19] Extracting length features."
#> [1] "[2024-11-17 15:13:20] Extracting motif features."
#> [1] "[2024-11-17 15:13:20] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-11-17 15:13:20] Extracting kmer features."
#> [1] "[2024-11-17 15:13:20] Extracting A Ratio features."
#> [1] "[2024-11-17 15:13:21] Saving Result"
#> [1] "[2024-11-17 15:13:22] Extracting A3SS features Finished"
#> [1] "[2024-11-17 15:13:22] Extracting A5SS features..."
#> [1] "[2024-11-17 15:13:22] Loading events..."
#> [1] "[2024-11-17 15:13:32] Parsing events region..."
#> [1] "[2024-11-17 15:13:40] Extracting length features."
#> [1] "[2024-11-17 15:13:40] Extracting motif features."
#> [1] "[2024-11-17 15:13:40] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-11-17 15:13:40] Extracting kmer features."
#> [1] "[2024-11-17 15:13:40] Extracting A Ratio features."
#> [1] "[2024-11-17 15:13:42] Saving Result"
#> [1] "[2024-11-17 15:13:42] Extracting A5SS features Finished"
#> [1] "[2024-11-17 15:13:42] Extracting MXE features..."
#> [1] "[2024-11-17 15:13:42] Loading events..."
#> [1] "[2024-11-17 15:13:52] Parsing events region..."
#> [1] "[2024-11-17 15:14:00] Extracting length features."
#> [1] "[2024-11-17 15:14:00] Extracting motif features."
#> [1] "[2024-11-17 15:14:00] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-11-17 15:14:01] Extracting kmer features."
#> [1] "[2024-11-17 15:14:01] Extracting A Ratio features."
#> [1] "[2024-11-17 15:14:01] Saving Result"
#> [1] "[2024-11-17 15:14:01] Extracting MXE features Finished"
#> [1] "[2024-11-17 15:14:01] Extracting RI features..."
#> [1] "[2024-11-17 15:14:01] Loading events..."
#> [1] "[2024-11-17 15:14:11] Parsing events region..."
#> [1] "[2024-11-17 15:14:18] Extracting length features."
#> [1] "[2024-11-17 15:14:18] Extracting motif features."
#> [1] "[2024-11-17 15:14:18] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-11-17 15:14:19] Extracting kmer features."
#> [1] "[2024-11-17 15:14:19] Extracting A Ratio features."
#> [1] "[2024-11-17 15:14:20] Saving Result"
#> [1] "[2024-11-17 15:14:20] Extracting RI features Finished"
#> [1] "[2024-11-17 15:14:20] Extracting SE features..."
#> [1] "[2024-11-17 15:14:20] Loading events..."
#> [1] "[2024-11-17 15:14:31] Parsing events region..."
#> [1] "[2024-11-17 15:14:42] Extracting length features."
#> [1] "[2024-11-17 15:14:44] Extracting motif features."
#> [1] "[2024-11-17 15:14:44] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2024-11-17 15:14:47] Extracting kmer features."
#> [1] "[2024-11-17 15:14:47] Extracting A Ratio features."
#> [1] "[2024-11-17 15:14:47] Saving Result"
#> [1] "[2024-11-17 15:14:47] Extracting SE features Finished"
#> [1] "[2024-11-17 15:14:47] step3 Combining events feature ======="
#> [1] "[2024-11-17 15:14:47] Parsing A3SS features..."
#> [1] "[2024-11-17 15:14:48] Parsing A3SS features Finished"
#> [1] "[2024-11-17 15:14:48] Parsing A5SS features..."
#> [1] "[2024-11-17 15:14:48] Parsing A5SS features Finished"
#> [1] "[2024-11-17 15:14:48] Parsing MXE features..."
#> [1] "[2024-11-17 15:14:48] Parsing MXE features Finished"
#> [1] "[2024-11-17 15:14:48] Parsing RI features..."
#> [1] "[2024-11-17 15:14:49] Parsing RI features Finished"
#> [1] "[2024-11-17 15:14:49] Parsing SE features..."
#> [1] "[2024-11-17 15:14:50] Parsing SE features Finished"
#> [1] "[2024-11-17 15:14:50] step4 Encoding events feature ======="
#> [1] "[2024-11-17 15:14:50] A3SS event encoding..."
#> 1/4 [======>.......................] - ETA: 0s 4/4 [==============================] - 0s 1ms/step
#> [1] "[2024-11-17 15:15:00] A3SS event encoding Finish."
#> [1] "[2024-11-17 15:15:00] A5SS event encoding..."
#> 1/4 [======>.......................] - ETA: 0s 4/4 [==============================] - 0s 1ms/step
#> [1] "[2024-11-17 15:15:08] A5SS event encoding Finish."
#> [1] "[2024-11-17 15:15:08] MXE event encoding..."
#> 1/2 [==============>...............] - ETA: 0s 2/2 [==============================] - 0s 2ms/step
#> [1] "[2024-11-17 15:15:16] MXE event encoding Finish."
#> [1] "[2024-11-17 15:15:16] RI event encoding..."
#>  1/12 [=>............................] - ETA: 0s 12/12 [==============================] - 0s 1ms/step
#> [1] "[2024-11-17 15:15:24] RI event encoding Finish."
#> [1] "[2024-11-17 15:15:24] SE event encoding..."
#>  1/61 [..............................] - ETA: 4s 50/61 [=======================>......] - ETA: 0s 61/61 [==============================] - 0s 1ms/step
#> [1] "[2024-11-17 15:15:34] SE event encoding Finish."
#> [1] "[2024-11-17 15:15:34] step5 Calculate splicing regulation distance and Combine distance ======="
#> [1] "384 rbps are used to calculate splicing regulation information"
#> [1] "Save data"
#> [1] "Save data Finished"
#> [1] "[2024-11-17 15:15:55] step6 Calculate combined event similarity ======="
#> [1] "[2024-11-17 15:15:57] Calculate  A3SS event Similarity"
#> 
#> Attaching package: 'Matrix'
#> The following object is masked from 'package:S4Vectors':
#> 
#>     expand
#> [1] "[2024-11-17 15:16:06] Calculate A3SS event Similarity Finished"
#> [1] "[2024-11-17 15:16:07] Calculate  A5SS event Similarity"
#> [1] "[2024-11-17 15:16:16] Calculate A5SS event Similarity Finished"
#> [1] "[2024-11-17 15:16:17] Calculate  MXE event Similarity"
#> [1] "[2024-11-17 15:16:26] Calculate MXE event Similarity Finished"
#> [1] "[2024-11-17 15:16:27] Calculate  RI event Similarity"
#> [1] "[2024-11-17 15:16:36] Calculate RI event Similarity Finished"
#> [1] "[2024-11-17 15:16:37] Calculate  SE event Similarity"
#> [1] "[2024-11-17 15:16:46] Calculate SE event Similarity Finished"
#> [1] "[2024-11-17 15:16:46] Calculate event similarity Finished."

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

A list of event similarity for different types of splicing events will
be saved to `work_path/imputation/event_similarity/event.similars.rds`;

A list of different types of cell similarity and a list of the number of
neighbors will be saved to rds file
to`work_path/imputation/cell_similarity/cell.similars.rds` and
`work_path/imputation/cell_similarity/dyk.cell.rds`

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

`alpha_cell`: random walk probability (1-restart probability).

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

`alpha_event`: random walk probability (1-restart probability)

`decay_event`: threshold of change in the similarity matrix

### Step6. Imputation

Based on these weighted similarity networks, SCSES next will use three
imputation strategies to aggregate the information across similar cells
or events to impute read count or PSI value

``` r
Imputed.data.path = ImputationAll(paras)
```

``` r
Imputed.data.path = ImputationAll(paras)
#> [1] "[2024-11-17 15:16:47] Get imputed result using cell similarity and event similarity."
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
#> [1] "[2024-11-17 15:16:48] Running Event_type=A3SS;cell_similarity_feature=PSI"
#> [1] "[2024-11-17 15:16:48] Save data"
#> [1] "[2024-11-17 15:16:48] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-499989986.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-499989986.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:17:01] Running Event_type=A3SS;cell_similarity_feature=RC"
#> [1] "[2024-11-17 15:17:01] Save data"
#> [1] "[2024-11-17 15:17:01] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-500010096.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-500010096.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:17:15] Running Event_type=A3SS;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-11-17 15:17:15] Save data"
#> [1] "[2024-11-17 15:17:15] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-499993818.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-499993818.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:17:28] Running Event_type=A5SS;cell_similarity_feature=PSI"
#> [1] "[2024-11-17 15:17:28] Save data"
#> [1] "[2024-11-17 15:17:28] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-499977711.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-499977711.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:17:41] Running Event_type=A5SS;cell_similarity_feature=RC"
#> [1] "[2024-11-17 15:17:41] Save data"
#> [1] "[2024-11-17 15:17:41] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-500001640.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-500001640.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:17:55] Running Event_type=A5SS;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-11-17 15:17:55] Save data"
#> [1] "[2024-11-17 15:17:55] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-500001365.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-500001365.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:18:08] Running Event_type=MXE;cell_similarity_feature=PSI"
#> [1] "[2024-11-17 15:18:08] Save data"
#> [1] "[2024-11-17 15:18:08] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-499999413.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-499999413.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:18:21] Running Event_type=MXE;cell_similarity_feature=RC"
#> [1] "[2024-11-17 15:18:21] Save data"
#> [1] "[2024-11-17 15:18:21] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-500012791.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-500012791.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:18:34] Running Event_type=MXE;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-11-17 15:18:34] Save data"
#> [1] "[2024-11-17 15:18:34] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-499980133.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-499980133.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:18:47] Running Event_type=RI;cell_similarity_feature=PSI"
#> [1] "[2024-11-17 15:18:47] Save data"
#> [1] "[2024-11-17 15:18:47] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-499992000.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-499992000.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:19:01] Running Event_type=RI;cell_similarity_feature=RC"
#> [1] "[2024-11-17 15:19:01] Save data"
#> [1] "[2024-11-17 15:19:01] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-500020303.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-500020303.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:19:14] Running Event_type=RI;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-11-17 15:19:14] Save data"
#> [1] "[2024-11-17 15:19:14] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-500012119.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-500012119.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:19:27] Running Event_type=SE;cell_similarity_feature=PSI"
#> [1] "[2024-11-17 15:19:27] Save data"
#> [1] "[2024-11-17 15:19:27] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-500008696.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-500008696.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:19:40] Running Event_type=SE;cell_similarity_feature=RC"
#> [1] "[2024-11-17 15:19:40] Save data"
#> [1] "[2024-11-17 15:19:40] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-500002877.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-500002877.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:19:53] Running Event_type=SE;cell_similarity_feature=EXP_RBP"
#> [1] "[2024-11-17 15:19:53] Save data"
#> [1] "[2024-11-17 15:19:53] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_data_2349247-499996445.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//imputation_result_2349247-499996445.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//mat_Imputation.log 2>&1"
#> [1] "[2024-11-17 15:20:06] Get imputed result using cell similarity and event similarity Finish."
Imputed_seperated = readRDS(Imputed.data.path)
str(Imputed_seperated,max.level=3)
#> List of 2
#>  $ cell      :List of 6
#>   ..$ PSI_PSI    : num [1:2550, 1:15] 0.881 0.955 0.211 0.138 0.721 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ PSI_RC     : num [1:2550, 1:15] 1 0.889 0.371 0.256 0.662 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ RC_PSI     : num [1:2550, 1:15] 0.8106 0.9394 0.261 0.0567 0.7259 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ RC_RC      : num [1:2550, 1:15] 1 0.849 0.422 0.177 0.611 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ EXP_RBP_PSI: num [1:2550, 1:15] 0.8759 0.9082 0.2418 0.0698 0.5931 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ EXP_RBP_RC : num [1:2550, 1:15] 1 0.844 0.414 0.192 0.61 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>  $ cell_event:List of 3
#>   ..$ PSI_PSI    : num [1:2550, 1:15] 0.566 0.684 0.486 0.553 0.697 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ RC_PSI     : num [1:2550, 1:15] 0.535 0.645 0.481 0.557 0.663 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ EXP_RBP_PSI: num [1:2550, 1:15] 0.536 0.601 0.502 0.554 0.629 ...
#>   .. ..- attr(*, "dimnames")=List of 2
# explain:
# For example:
# cell/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw PSI with cell similarities.
# cell/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw inclusion and exclusion read counts with cell similarities, and then calculate the imputed PSI.
# cell_event/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw inclusion and exclusion read counts with cell similarities, and then calculate the imputed PSI. Further impute the results using event similarities.
```

Results of each imputation strategy will be saved
into`work_path/imputation`.

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
```

``` r
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
Imputed.data.final.path = Estimation(paras,rds_imputed_file = Imputed.data.path)
#> [1] "[2024-11-17 15:20:07] Combine imputed psi."
#> [1] "[2024-11-17 15:20:07] Loading data..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_seperated_500012042.rds"
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
#> [1] "[2024-11-17 15:20:10] Combine imputed psi Finish."
Imputed_combined = readRDS(Imputed.data.final.path)
str(Imputed_combined,max.level=2)
#> List of 3
#>  $ PSI    : num [1:2550, 1:15] 1 0.955 0.211 0.138 0.721 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ RC     : num [1:2550, 1:15] 1 0.9394 0.261 0.0567 0.7259 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ EXP_RBP: num [1:2550, 1:15] 1 0.9082 0.2418 0.0698 0.5931 ...
#>   ..- attr(*, "dimnames")=List of 2
# The finally imputed PSI were named by cell-cell splicing similarity features, including raw event PSI(PSI), and raw junction read counts(RC), and RBP expression(EXP_RBP).
```

A list of final imputation of PSI values will be saved by serialized R
object format (.rds) in `work_path/imputation/Imputed_combined*`, where
**\*** is a random number representing different execution. The `.rds`
file can be loaded in R environment by `readRDS` function.

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
ftpsi.path = getFtRawPSI(paras)
ftrds.path = mergeFtSplicingValue(paras)
ftmodel.path = FtClassifier(paras)
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
ImputedFt.data.final.path = Estimation(paras,rds_imputed_file = Imputed.data.path)
```

``` r
ftrc.path = getFtRawRC(paras)
#> [1] "Loading splicing events for classifer fine tune..."
#> [1] "Checking cells..."
#> [1] "15 cells are considered."
#> [1] "Output: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/splicing_value_ft/"
#> [1] "[2024-11-17 15:20:10] Counting reads of A3SS events..."
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
#> [1] "[2024-11-17 15:20:13] Counting reads of A3SS events Finish."
#> [1] "[2024-11-17 15:20:13] Counting reads of A5SS events..."
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
#> [1] "[2024-11-17 15:20:15] Counting reads of A5SS events Finish."
#> [1] "[2024-11-17 15:20:15] Counting reads of MXE events..."
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
#> [1] "[2024-11-17 15:20:16] Counting reads of MXE events Finish."
#> [1] "[2024-11-17 15:20:16] Counting reads of SE events..."
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
#> [1] "[2024-11-17 15:20:19] Counting reads of SE events Finish."
ftpsi.path = getFtRawPSI(paras)
#> [1] "Checking raw reads..."
#> [1] "Loading splicing events for classifer fine tune..."
#> [1] "[2024-11-17 15:20:19] Calculating PSI value of A3SS events..."
#> [1] "[2024-11-17 15:20:19] Calculating PSI value of A3SS events Finish."
#> [1] "[2024-11-17 15:20:19] Calculating PSI value of A5SS events..."
#> [1] "[2024-11-17 15:20:19] Calculating PSI value of A5SS events Finish."
#> [1] "[2024-11-17 15:20:19] Calculating PSI value of MXE events..."
#> [1] "[2024-11-17 15:20:19] Calculating PSI value of MXE events Finish."
#> [1] "[2024-11-17 15:20:19] Calculating PSI value of SE events..."
#> [1] "[2024-11-17 15:20:19] Calculating PSI value of SE events Finish."
ftrds.path = mergeFtSplicingValue(paras)
ftmodel.path = FtClassifier(paras)
#> [1] "Reading true Ft PSI..."
#> [1] "Loading Pre-training classifer..."
#> [1] "[2024-11-17 15:20:19] Classifer fine tune"
#> [1] "[2024-11-17 15:20:19] Processing raw Ft data..."
#> [1] "Checking data..."
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=PSI;RC;EXP_RBP  checked"
#> [1] "Calculating Classifier Features..."
#> [1] "[2024-11-17 15:20:21] Save data"
#> [1] "[2024-11-17 15:20:21] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499990180.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499990180.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:20:44] Save data"
#> [1] "[2024-11-17 15:20:44] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499962078.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499962078.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:21:07] Save data"
#> [1] "[2024-11-17 15:21:07] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499969972.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499969972.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:21:30] Save data"
#> [1] "[2024-11-17 15:21:30] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499988493.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499988493.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:21:53] Save data"
#> [1] "[2024-11-17 15:21:53] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499977550.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499977550.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:22:16] Save data"
#> [1] "[2024-11-17 15:22:16] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499987047.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499987047.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:22:39] Save data"
#> [1] "[2024-11-17 15:22:39] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499984279.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499984279.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:23:03] Save data"
#> [1] "[2024-11-17 15:23:03] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-500001823.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-500001823.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:23:28] Save data"
#> [1] "[2024-11-17 15:23:28] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499980114.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499980114.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:23:53] Save data"
#> [1] "[2024-11-17 15:23:53] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-500005357.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-500005357.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:24:16] Save data"
#> [1] "[2024-11-17 15:24:16] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-500014346.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-500014346.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:24:40] Save data"
#> [1] "[2024-11-17 15:24:40] Save data Finished"
#> [1] "bash /tmp/RtmppACzdc/temp_libpath1e705f31c80257/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_data_2349247-499967212.h5 /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//imputation1_result_2349247-499967212.mat >> /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2024-11-17 15:25:03]  Model training;similarity_type=PSI"
#> [1] "[2024-11-17 15:25:04]  Model training;similarity_type=RC"
#> [1] "[2024-11-17 15:25:05]  Model training;similarity_type=EXP_RBP"
#> [1] "[2024-11-17 15:25:05] Classifer fine tune Finish."
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
ImputedFt.data.final.path = Estimation(paras,rds_imputed_file = Imputed.data.path)
#> [1] "[2024-11-17 15:25:05] Combine imputed psi."
#> [1] "[2024-11-17 15:25:05] Loading data..."
#> [1] "Input: /disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_seperated_500012042.rds"
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
#> [1] "[2024-11-17 15:25:09] Combine imputed psi Finish."

print(ImputedFt.data.final.path)
#> [1] "/disk/lvxuan/Single-Splicing/result/cell_line/scses_test/imputation//Imputed_combined_499987971.rds"
```

The final imputation results can be loaded by:

``` r
Imputed_combined = readRDS(ImputedFt.data.final.path)
Imputed_combined[["EXP_RBP"]][1:3,1:3]
#>                                                                                                                                                                                                                                        SRR11826368.bam
#> isoform1=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178086787:+@exon:chr2:178086788-178088686:+|isoform2=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178084135:+@exon:chr2:178084136-178088686:+|HNRNPA3|A3SS       1.0000000
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS         0.9081793
#> isoform1=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041195:+@exon:chr3:184041196-184041381:+|isoform2=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041192:+@exon:chr3:184041193-184041381:+|EIF4G1|A3SS        0.2417979
#>                                                                                                                                                                                                                                        SRR11826371.bam
#> isoform1=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178086787:+@exon:chr2:178086788-178088686:+|isoform2=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178084135:+@exon:chr2:178084136-178088686:+|HNRNPA3|A3SS       0.9325028
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS         0.8759136
#> isoform1=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041195:+@exon:chr3:184041196-184041381:+|isoform2=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041192:+@exon:chr3:184041193-184041381:+|EIF4G1|A3SS        0.1922552
#>                                                                                                                                                                                                                                        SRR11826409.bam
#> isoform1=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178086787:+@exon:chr2:178086788-178088686:+|isoform2=exon:chr2:178083975-178084042:+@junction:chr2:178084043-178084135:+@exon:chr2:178084136-178088686:+|HNRNPA3|A3SS       0.9437169
#> isoform1=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025540:+@exon:chr4:146025541-146025667:+|isoform2=exon:chr4:146019084-146019571:+@junction:chr4:146019572-146025537:+@exon:chr4:146025538-146025667:+|ABCE1|A3SS         0.8796700
#> isoform1=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041195:+@exon:chr3:184041196-184041381:+|isoform2=exon:chr3:184040871-184041029:+@junction:chr3:184041030-184041192:+@exon:chr3:184041193-184041381:+|EIF4G1|A3SS        0.2184908
```

The reference results of the test data can be found
[here](https://github.com/lvxuan12/SCSES/blob/main/analysis/Imputed_combined_499993907.rds),
or loaded by command:

``` r
ref.path <- paste0(system.file("analysis", package = "SCSES"),'/Imputed_combined_499993907.rds')
ref.result <- readRDS(ref.path)
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

annotation=read.table(paste0(system.file("analysis", package = "SCSES"),'/cell_line_annotation.txt'),sep="\t")
# annotation=read.table("/disk/lvxuan/Single-Splicing/src/package/SCSES/analysis/cell_line_annotation.txt",sep="\t")
Imputed_combined=readRDS(paste0(system.file("analysis", package = "SCSES"),'/Imputed_combined_499993907.rds'))
# Imputed_combined=readRDS("/disk/lvxuan/Single-Splicing/src/package/SCSES/analysis/Imputed_combined_499993907.rds")
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

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />
