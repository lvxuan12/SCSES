Menu
================

- [SCSES](#scses)
  - [Hardware requirements](#hardware-requirements)
  - [Software requirements](#software-requirements)
    - [OS Requirements](#os-requirements)
  - [Installation](#installation)
    - [Installation with docker file](#installation-with-docker-file)
      - [1. Install Docker Client:](#1-install-docker-client)
      - [2. Download Dockerfile](#2-download-dockerfile)
      - [3. Build Docker Image](#3-build-docker-image)
      - [4. Create Docker Container](#4-create-docker-container)
      - [5. Access RStudio Server](#5-access-rstudio-server)
    - [Installation with conda](#installation-with-conda)
      - [Step 1: Environment Setup](#step-1-environment-setup)
      - [Step 2: Install Python Dependencies](#step-2-install-python-dependencies)
      - [Step 3: Install Tools](#step-3-install-tools)
      - [Step 4: Install SCSES R Package](#step-4-install-scses-r-package)
      - [Tips for some Installation errors](#tips-for-some-installation-errors)
  - [SCSES input](#scses-input)
    - [1. BAM Files](#1-bam-files)
    - [2. Reference Genome Files](#2-reference-genome-files)
    - [3. Configuration File](#3-configuration-file)
    - [4. Phast conservation file in bigWig format](#4-phast-conservation-file-in-bigwig-format)
    - [5. RBP](#5-rbp)
  - [Getting started](#getting-started)
    - [Download Test Data](#download-test-data)
      - [Setup](#setup)
    - [Step-by-Step Analysis](#step-by-step-analysis)
      - [Step 1. Read config file](#step-1-read-config-file)
      - [Step 2. Get gene expression](#step-2-get-gene-expression)
      - [Step 3. Detect splicing events](#step-3-detect-splicing-events)
      - [Step 4. Quantify splicing events](#step-4-quantify-splicing-events)
      - [Step 5. Constructs similarity networks](#step-5-constructs-similarity-networks)
      - [Step 6. Imputation](#step-6-imputation)
      - [Step 7. Estimation](#step-7-estimation)
      - [Step 8. Cell Clustering](#step-8-cell-clustering)

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

This package is supported for Linux.

The package has been tested on Debian-11.21.

## Installation

### Installation with docker file

SCSES provides a Docker-based installation method to simplify the setup
of all dependencies and requirements. Please follow the steps below to
build the Docker image and start a container to use SCSES:

#### 1. Install Docker Client:

Please install the [Docker
client](https://www.docker.com/products/docker-desktop) on the host
machine.

#### 2. Download Dockerfile

The Dockerfile of SCSES can be downloaded from:
<https://github.com/lvxuan12/SCSES/blob/main/SCSES.dockerfile>.

#### 3. Build Docker Image

You can build the SCSES Docker image using the command:

``` bash
docker build -t scses -f SCSES.dockerfile .
```

#### 4. Create Docker Container

After building the image, create a Docker container with the following
command:

``` bash
docker run -d -p [exported port]:8787 -e PASSWORD=[user password] -v [local directory]:/data --name test scses
```

`[exported port]`: An **unused** port on the host machine to access the container. To show all **used** ports on the host machine, 
input the following command in linux terminal `netstat -tuwanp 2|awk '{print $4}'|cut -d ":" -f 2|sort|uniq -c` 
or ` Get-NetTCPConnection | Where-Object { $_.State -eq "Listen" } | Select-Object -ExpandProperty LocalPort|Sort-Object | Group-Object | Select-Object -Property Count, Name` in Windows Powershell.

`[user password]`: A user-defined password for logging into the RStudio
server.

`[local directory]`: A local directory mapped to the container for data
storage and sharing.

#### 5. Access RStudio Server

Now, you can access the RStudio server by opening a web browser and
navigating to `[host IP]:[exported port]`. The username to log in Rstudio server is
**`rstudio`** and the password is use-defined in the `docker run` command.

In this pre-configd RStudio server environment, SCSES and all its
dependencies are correctly installed and ready for use.

Please refer to [Getting started](#getting-started) to start the first
experience with SCSES.

Enjoy!

### Installation with conda

#### Step 1: Environment Setup

We recommend a new **conda** environment to install SCSES:

``` bash
conda create -n SCSES_test python=3.11
conda activate SCSES_test
```

To use SCSES, you will need to install R, Python, Matlab Compiler
Runtime(v9.13), and Java(v17.0.10).

``` bash
## install R in conda environment
conda install -c conda-forge r-base=4.3.1
```

The MCR is quite large, so downloading may take some time.

``` bash
## install MCR
mkdir /path/to/MCR && \
cd /path/to/MCR && \
wget https://ssd.mathworks.com/supportfiles/downloads/R2022b/Release/10/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2022b_Update_10_glnxa64.zip && \
unzip -q MATLAB_Runtime_R2022b_Update_10_glnxa64.zip && \
./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent
```

#### Step 2: Install Python Dependencies

``` bash
pip install pandas numpy scipy scikit-learn
pip install keras==2.15.0
pip install tensorflow==2.15.0.post1
```

#### Step 3: Install Tools

To detect splicing events You will need to install rMATS, MAJIQ,
IRFinder. rMATS should be built in the same environment with SCSES (same
python). MAJIQ should be built in a new environment due to conflict of
python package version，

##### 3.1 [rMATS](https://github.com/Xinglab/rmats-turbo)

``` bash
wget https://github.com/Xinglab/rmats-turbo/releases/download/v4.3.0/rmats_turbo_v4_3_0.tar.gz
tar -zxvf rmats_turbo_v4_3_0.tar.gz
cd rmats_turbo_v4_3_0/
# Install dependencies
pip install Cython
# Add to PATH
./build_rmats
export PATH=/path/to/rmats_turbo_v4_3_0/:$PATH
```

##### 3.2 [MAJIQ](https://biociphers.bitbucket.io/majiq-docs/index.html)

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

##### 3.3 [IRFinder](https://github.com/dgaolab/IRFinder)

``` bash
wget https://github.com/RitchieLabIGH/IRFinder/archive/refs/tags/v2.0.1.tar.gz
tar -zxvf v2.0.1.tar.gz
export PATH=/path/to/IRFinder-2.0.1/bin/:$PATH
```

**NOTE:** [STAR](https://github.com/alexdobin/STAR) is required to build
IRFinder reference! To run IRFinder correctly, you also need to install
STAR.

##### 3.4 [samtools](https://github.com/samtools/samtools)

#### Step 4: Install SCSES R Package

Currently SCSES can only be installed from GitHub.

To install SCSES, type the following command in **R**:

``` r
install.packages("remotes")
remotes::install_version("Matrix", version = "1.6-5")
install.packages("curl",config.vars='LIB_DIR=/usr/lib/x86_64-linux-gnu/pkgconfig/')
options(download.file.method = "wget", times=100)
remotes::install_github("lvxuan12/SCSES")
```

#### Tips for some Installation errors

##### Issue 1. cannot find fftw.h

``` bash
conda install conda-forge::fftw
# config FFTW
export FFTW_CFLAGS=" -I/path/to/miniconda3/envs/scses/include/"
export FFTW_LIBS=" -L/path/to/miniconda3/envs/scses/lib -lfftw3"
```

##### Issue 2. cannot find -lxml2

``` bash
conda install conda-forge::libxml2
```

##### Issue 3. cannot find -lsz

``` bash
ln -s /usr/lib/x86_64-linux-gnu/libsz.so /path/to/miniconda/envs/SCSES_test/lib/libsz.so
```

##### Issue 4. error: ‘::timespec_get’ has not been declared

``` bash
conda upgrade -c conda-forge --all 
```

## SCSES input

SCSES requires five essential input files:

### 1. BAM Files

- **Format**: Coordinate-sorted BAM files with index (.bai)

### 2. Reference Genome Files

- **FASTA**: Reference genome sequence
- **GTF**: Gene annotations
- **GFF3**: Gene annotations

### 3. Configuration File

SCSES requires a json-based configuration file to set all parameters in the algorithm. Here is a [demo](https://github.com/lvxuan12/SCSES/blob/main/inst/analysis/cell_line.json) of the configure file.
For a detailed explanation of the configuration file, please refer to the [ConfigurationGuide.txt](https://github.com/lvxuan12/SCSES/blob/main/ConfigurationGuide.txt).

SCSES provides a shiny app to help you to generate the confugre file. You can start the app by `createConfigshiny` function.

If you use it in the SCSES docker container, the command should be :

``` r
library(SCSES)
createConfigshiny(host = "localhost",launch.browser=TRUE) 
```

For non-docker users, the full command should be:
``` r
library(SCSES)
createConfigshiny(host, port, launch.browser=FALSE) 
```
You should set the following parameters:

- host: the server's IP address, for local access, you can set as "localhost" or "127.0.0.1"

- port: The TCP port that the application should listen on. If the port is not specified, and the shiny.port option is set (with options(shiny.port = XX)), then that port will be used. Otherwise, use a random port between 3000:8000, excluding ports that are blocked by Google Chrome for being considered unsafe: 3659, 4045, 5060, 5061, 6000, 6566, 6665:6669 and 6697. Up to twenty random ports will be tried.

- launch.borwser: if launch the app in the default web browser automatically, default is FALSE. Setting launch.browser = TRUE may cause errors in headless environments (servers without GUI) or when no default browser is configured

After running `createConfigshiny`, you will see a URL appear in
the console. Copy this URL and paste it into your web browser to access
the application.

An interactive window will popup, which allow you to
fill some parameters, such as Bam File Path, and Work Path. The meaning of each parameter can be found by hovering the mouse over the widget.

Finally, you can click `Create Config` button and a json file will be generated in the `work_path` you provided if successful.

If you are using **test data** in this Tutorial, you should use `createDemoConfigshiny`
function instead to build the configuration file:

``` r
library(SCSES)
createDemoConfigshiny(host = "localhost", launch.browser=TRUE) 
```

**Note**: The test dataset includes fewer cells and chromosomes to ensure faster completation of the Tutorial. 
Therefore, the default parameters in `createConfigshiny` are not suitable. Please use the `createDemoConfigshiny` function instead, 
which provides default values optimized for the test dataset.


### 4. Phast conservation file in bigWig format

For human and mouse, you could download it directly from UCSC browser:

| Species | File | Size | Link |
|----|----|----|----|
| Human (hg38) | hg38.phastCons100way.bw | 5.5 GB | [Download](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw) |
| Human (hg19) | hg19.100way.phastCons.bw | 5.4 GB | [Download](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw) |
| Mouse (mm10) | mm10.60way.phastCons.bw | 4.3 GB | [Download](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw) |

### 5. RBP

Genes annotated as RBP are required to constructs similarity networks.
For human and mouse, you could read from package:

``` r
# Load human RBP list
rbp_human <- system.file("extdata/rbp/human_rbp.txt", package = "SCSES")
rbp_list <- readLines(rbp_human)
cat("Total RBPs:", length(rbp_list), "\n")
#> Total RBPs: 1456
cat("First 10 RBPs:", head(rbp_list, 10), sep = "\n")
#> First 10 RBPs:
#> A1CF
#> AC004381.6
#> ACIN1
#> ACO1
#> AKAP1
#> ALKBH8
#> ALYREF
#> ANKHD1
#> ANKRD17
#> APTX
# Load mouse RBP list
rbp_mouse <- system.file("extdata/rbp/mouse_rbp.txt", package = "SCSES")
rbp_list <- readLines(rbp_mouse)
cat("Total RBPs:", length(rbp_list), "\n")
#> Total RBPs: 611
cat("First 10 RBPs:", head(rbp_list, 10), sep = "\n")
#> First 10 RBPs:
#> Mcts1
#> Mkrn2
#> Hnrnpd
#> Fmr1
#> Snrpn
#> Pcbp3
#> Trmt1
#> Supt6h
#> Cugbp2
#> Sf3a1
```

## Getting started

### Download Test Data

This dataset includes BAM files for three cell lines (HCT116, HepG2,
HL-60), each containing five cells and other input files that are
essential for running the SCSES package.

**Download**: <https://doi.org/10.5281/zenodo.15688700>

| File Type | File Name | Description |
|----|----|----|
| **BAM** | `*.bam` | BAM files for three cell lines, each containing five cells |
| **BAI** | `*.bam.bai` | BAM index files |
| **TXT** | `annotation.txt` | Cell identities |
| **FASTA** | `test.fa` | Reference genome sequence |
| **FAI** | `test.fa.fai` | Reference genome sequence index |
| **Annotation** | `test.gtf` | Gene annotation file |
| **Annotation** | `test.gff3` | Gene annotation file |
| **TXT** | `human_rbp.txt` | RNA-binding proteins list |
| **PhastCons** | `test_phastCons.bw` | Conservation scores |
| **JSON** | `cell_line.json` | Parameters config file |

#### Setup

After download, ensure you have:

- 15 BAM files + index files

- annotation.txt: cell identities

- test.fa and test.fai: reference genome sequence

- test.gtf and test.gff: gene annotation file

- test_phastCons.bw: conservation scores

- human_rbp.txt: RBP genes list

- cell_line.json: parameters config file

Move 15 BAM files and their index to `bam` directory.

Move other input data to `refgenome` directory.

``` bash
# Example:
ls /disk/share/lvxuan/SCSES_test/bam/
ls /disk/share/lvxuan/SCSES_test/refgenome/
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
#> annotation.txt
#> cell_line.json
#> human_rbp.txt
#> test.fa
#> test.fa.fai
#> test.gff3
#> test.gtf
#> test_phastCons.bw
```

### Step-by-Step Analysis

#### Step 0. Get the cofigure file
The paramteter configuration for test dataset is the *cell_line.json* in the downloaded files. 

Alternatively, you can create the configuration file using the Shiny app.
**Note**: The test dataset includes fewer cells and chromosomes to ensure faster completation of the Tutorial. 
Therefore, the default parameters in `createConfigshiny` are not suitable. Please use the `createDemoConfigshiny` function instead, 
which provides default values optimized for the test dataset.

#### Step 1. Read config file

The `cell_line.json` file was downloaded previously

``` r
## Loading packages
library(SCSES)
#> Loading required package: BiocManager
#> Bioconductor version '3.18' is out-of-date; the current release version '3.21'
#>   is available with R version '4.5'; see https://bioconductor.org/install
#> Warning: replacing previous import 'hdf5r::h5version' by 'rhdf5::h5version'
#> when loading 'SCSES'
#> Warning: replacing previous import 'hdf5r::h5const' by 'rhdf5::h5const' when
#> loading 'SCSES'
#> Warning: replacing previous import 'fs::path' by 'rtracklayer::path' when
#> loading 'SCSES'
#> Warning: replacing previous import 'hdf5r::values' by 'rtracklayer::values'
#> when loading 'SCSES'

# Load configuration file
config_file <- system.file("analysis/cell_line.json", package = "SCSES")
paras <- readSCSESconfig(config_file)

# Verify configuration
cat("Dataset:", paras$DataSet, "\n")
#> Dataset: cell_line
cat("BAM path:", paras$Basic$bam_path, "\n")
#> BAM path: /disk/share/lvxuan/SCSES_test/bam/
# all outputs will be saved in work_path
cat("Work path:", paras$Basic$work_path, "\n")
#> Work path: /disk/share/lvxuan/SCSES_test/
```

#### Step 2. Get gene expression

##### TPM matrix (for smart-seq2 dataset)

The TPM matrix of gene expression can be obtained by different methods.
We used
[featureCounts](https://github.com/ShiLab-Bioinformatics/subread) to
generate this matrix.

You can use `getGeneExpression` to run `featureCounts`, which will save
featureCounts output to `work_path/expr/` and `getEXPmatrix` to generate
TPM matrix, which will save gene expression count and TPM matrix to
`work_path/rds/`.

``` r
cat("Quantifying gene expression...\n")
#> Quantifying gene expression...
featurecounts.path = getGeneExpression(paras) 
#> [1] "[2025-06-21 10:33:07.476851] Detect gene expression: bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/shell/run_featurecounts.sh /disk/share/lvxuan/SCSES_test//expr/ /disk/share/lvxuan/SCSES_test/refgenome/test.fa /disk/share/lvxuan/SCSES_test/refgenome/test.gtf /disk/share/lvxuan/SCSES_test/bam/ 20 cell_line paired /disk/software/subread-2.0.6-source/bin/featureCounts >> /disk/share/lvxuan/SCSES_test//runfeatureCounts.log 2>&1"
#> [1] "[2025-06-21 10:33:18.165756] Detect gene expression Finish."
rds.path = getEXPmatrix(paras)

# Load and examine results
print(rds.path)
#> [1] "/disk/share/lvxuan/SCSES_test//rds/"
list.files(rds.path)
#> [1] "count_norm.rds" "count.rds"      "event.rds"      "psi.rds"       
#> [5] "rc.rds"
tpm = readRDS(paste0(rds.path,'/count_norm.rds'))
tpm[1:5,1:5]
#>            SRR11826368.bam SRR11826371.bam SRR11826409.bam SRR11826436.bam
#> DDX11L1           0.000000        0.000000        0.000000               0
#> WASH7P            3.402929        3.815272        5.394269               0
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

##### Normalized UMI count matrix (for UMI based dataset)

You can use `get10XEXPmatrix` to generate Normalized UMI count matrix
from 10X CellRanger hdf5 file, which will save normalized UMI count to
`work_path/rds/`.

``` r
# Seurat package is required

rds.path <- get10XEXPmatrix(paras, expr_path)
```

#### Step 3. Detect splicing events

To define a global set of all splicing events, SCSES firstly merges all
bam files from every single cell to construct a pseudo-bulk bam file,
and identifies all types of splicing events by conventional algorithms.

##### For Smart-seq2 dataset

``` r
# Create pseudobulk for event detection

pseudobulk.path = createPseudobulk(paras)
#> [1] "Input: /disk/share/lvxuan/SCSES_test/bam/"
#> [1] "Output: /disk/share/lvxuan/SCSES_test//data/"
#> [1] "[2025-06-21 10:33:18.522907] Creating Pseudobulk directory..."
#> Warning in dir.create(path = pseudobulk.path, recursive = T):
#> '/disk/share/lvxuan/SCSES_test//data' already exists
#> [1] "Pseudobulk bam file exists."
#> [1] "[2025-06-21 10:33:18.526561] Merge Bam Files: /disk/software/samtools/bin/samtools merge -f -@ 20 /disk/share/lvxuan/SCSES_test//data//all.bam /disk/share/lvxuan/SCSES_test/bam//*.bam --no-PG"
#> [1] "[2025-06-21 10:33:22.380382] Merge Bam Files Finish."
#> [1] "[2025-06-21 10:33:22.381074] Bam File Index: /disk/software/samtools/bin/samtools index -@ 20 /disk/share/lvxuan/SCSES_test//data//all.bam"
#> [1] "[2025-06-21 10:33:23.04892] Bam File Index Finish."
print(pseudobulk.path)
#> [1] "/disk/share/lvxuan/SCSES_test//data/"
list.files(pseudobulk.path)
#> [1] "all.bam"     "all.bam.bai"

# Detect all AS events
#if you meet:
#irfinder: error while loading shared libraries: libboost_iostreams.so.1.71.0: cannot open shared object file: No such file or directory
#libboost_iostreams.so.1.71.0 exists in /disk/lvxuan/lib
old.ld=Sys.getenv("LD_LIBRARY_PATH")
Sys.setenv(LD_LIBRARY_PATH = paste0("/disk/lvxuan/lib:", old.ld))

event.path = detectEvents(paras)
#> [1] "Checking cells..."
#> [1] "Checking events..."
#> [1] "event_type=SE;RI;A3SS;A5SS;MXE  checked"
#> [1] "Output: /disk/share/lvxuan/SCSES_test//events/"
#> [1] "[2025-06-21 10:33:23.060435] Creating events directory..."
#> [1] "[2025-06-21 10:33:23.063627] MAJIQ..."
#> [1] "[2025-06-21 10:33:23.065012] Run MAJIQ: bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/shell/run_majiq.sh /disk/share/lvxuan/SCSES_test//events//majiq/ /disk/share/lvxuan/SCSES_test//data/ /disk/share/lvxuan/SCSES_test/refgenome/test.gff3 20 hg19 25 MAJIQ /disk/lvxuan/software/miniconda/bin /home/Liulab/wenx/majiq_license_academic_official.lic >> /disk/share/lvxuan/SCSES_test//events//detectEvents.log 2>&1"
#> [1] "[2025-06-21 10:34:22.754712] Run MAJIQ Finish."
#> [1] "[2025-06-21 10:34:22.756686] IRFinder..."
#> [1] "[2025-06-21 10:34:22.761278] Run IRFinder: bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/shell/run_irfinder.sh /disk/share/lvxuan/SCSES_test//events//IRFinder/ /disk/share/lvxuan/SCSES_test/refgenome/test.fa /disk/share/lvxuan/SCSES_test/refgenome/test.gtf /disk/share/lvxuan/SCSES_test//data/ 20 100 /disk/lvxuan/software/IRFinder-2.0.1/bin/IRFinder /disk/software/samtools/bin/samtools /usr/local/bin/STAR  >> /disk/share/lvxuan/SCSES_test//events//detectEvents.log 2>&1"
#> [1] "/disk/share/lvxuan/SCSES_test//events//IRFinder//all/IRFinder-IR-nondir.txt"
#> [1] "[2025-06-21 10:42:43.185548] Run IRFinder Finish."
#> [1] "[2025-06-21 10:42:43.186025] rMats..."
#> [1] "[2025-06-21 10:42:43.188481] Run rMats: bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/shell/run_rmats.sh /disk/share/lvxuan/SCSES_test//events//rMats/ /disk/share/lvxuan/SCSES_test//data/ /disk/share/lvxuan/SCSES_test/refgenome/test.gtf paired 100 20 /disk/lvxuan/software/SCSES_test/rmats_turbo_v4_3_0/rmats.py /disk/lvxuan/software/miniconda/envs/SCSES_test/bin/python3.11 >> /disk/share/lvxuan/SCSES_test//events//detectEvents.log 2>&1"
#> [1] "[2025-06-21 10:43:07.965503] Run rMats Finish."
#> [1] "[2025-06-21 10:43:07.967128] Generating SE event id"
#> arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
#> [1] "[2025-06-21 10:43:08.183085] Generating SE event id Finish."
#> [1] "[2025-06-21 10:43:08.184156] Generating RI event id"
#> arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
#> [1] "[2025-06-21 10:43:14.956387] Generating RI event id Finish."
#> [1] "[2025-06-21 10:43:14.957422] Generating A3SS event id"
#> [1] "[2025-06-21 10:43:15.095208] Generating A3SS event id Finish."
#> [1] "[2025-06-21 10:43:15.096094] Generating A5SS event id"
#> [1] "[2025-06-21 10:43:15.20304] Generating A5SS event id Finish."
#> [1] "[2025-06-21 10:43:15.203696] Generating MXE event id"
#> [1] "[2025-06-21 10:43:15.224238] Generating MXE event id Finish."
#> [1] "Total splicing event: SE=1685 RI=187 A3SS=88 A5SS=69 MXE=104"

# Summary of detected events
print(event.path)
#> [1] "/disk/share/lvxuan/SCSES_test//events/"
event_files <- list.files(event.path, pattern = "*.txt")
for(file in event_files) {
  events <- readLines(file.path(event.path, file))
  cat(gsub(".txt", "", file), "events:", length(events), "\n")
}
#> A3SS events: 85 
#> A5SS events: 69 
#> MXE events: 104 
#> RI events: 187 
#> SE events: 1685
```

Different types of splicing events will be saved to `work_path/events/`,
separately.

##### for UMI dataset

SCSES requires single cell bam files being saved in a directory. For
UMI-based dataset using CellRanger for data process, the function
`split10XBAM` can be used to get single cell bam files for each sample.

``` r
# CellRanger_path: directory to CellRanger output
# out_path: directory to save single cell bam
# java_path: directory to java
# core: the number of threads

splitbam.path = split10XBAM(CellRanger_path,out_path,java_path,core)

# path to single cell bam files should be added to bam_path in config file
paras$Basic$bam_path = splitbam.path
# pseudobulk.path = createPseudobulk(paras)
# It is not necessary to execute `createPseudobulk`, and the `possorted_genome_bam.bam`,`possorted_genome_bam.bam.bai` from `CellRanger_path` can be moved to `work_path/data/all.bam`.
event.path = detectEvents(paras)
```

#### Step 4. Quantify splicing events

According to splicing events detected in the previous step, SCSES then
quantify raw reads associated with these splicing events in each cell,
and construct the raw read count matrix and calculate the raw PSI
matrix.

The definition of PSI of different AS events: ![The definition of PSI of
different AS events:](png/PSI.png)

``` r
# Raw read counts and PSI calculation
rawrc.path = getRawRC(paras)
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "Checking cells..."
#> [1] "15 cells are considered."
#> [1] "Output: /disk/share/lvxuan/SCSES_test//splicing_value/"
#> [1] "Splicing event types: A3SS A5SS MXE RI SE"
#> [1] "[2025-06-21 10:43:15.26883] Counting reads of A3SS events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value//A3SS_rjm' already exists
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
#> [1] "[2025-06-21 10:43:55.124573] Counting reads of A3SS events Finish."
#> [1] "[2025-06-21 10:43:55.124884] Counting reads of A5SS events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value//A5SS_rjm' already exists
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
#> [1] "[2025-06-21 10:44:19.834326] Counting reads of A5SS events Finish."
#> [1] "[2025-06-21 10:44:19.834663] Counting reads of MXE events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value//MXE_rjm' already exists
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
#> [1] "[2025-06-21 10:44:52.104853] Counting reads of MXE events Finish."
#> [1] "[2025-06-21 10:44:52.105268] Counting reads of RI events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value//RI_rjm' already exists
#> [1] "Reading RJC File Progress: 0%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 10%"
#> [1] "Reading RJC File Progress: 20%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 30%"
#> [1] "Reading RJC File Progress: 40%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 50%"
#> [1] "Reading RJC File Progress: 60%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 70%"
#> [1] "Reading RJC File Progress: 80%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress: 90%"
#> [1] "Reading RJC File Progress:100%"
#> [1] "[2025-06-21 10:45:26.836598] Counting reads of RI events Finish."
#> [1] "[2025-06-21 10:45:26.837156] Counting reads of SE events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value//SE_rjm' already exists
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
#> [1] "[2025-06-21 10:46:12.432916] Counting reads of SE events Finish."
rawpsi.path = getRawPSI(paras)
#> [1] "Checking raw reads..."
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "[2025-06-21 10:46:12.4353] Calculating PSI value of A3SS events..."
#> [1] "[2025-06-21 10:46:12.439237] Calculating PSI value of A3SS events Finish."
#> [1] "[2025-06-21 10:46:12.439812] Calculating PSI value of A5SS events..."
#> [1] "[2025-06-21 10:46:12.442417] Calculating PSI value of A5SS events Finish."
#> [1] "[2025-06-21 10:46:12.44266] Calculating PSI value of MXE events..."
#> [1] "[2025-06-21 10:46:12.448437] Calculating PSI value of MXE events Finish."
#> [1] "[2025-06-21 10:46:12.448707] Calculating PSI value of RI events..."
#> [1] "[2025-06-21 10:46:12.453605] Calculating PSI value of RI events Finish."
#> [1] "[2025-06-21 10:46:12.453859] Calculating PSI value of SE events..."
#> [1] "[2025-06-21 10:46:12.500773] Calculating PSI value of SE events Finish."
# Merge and preprocess data
rawrds.path = mergeSplicingValue(paras)
processed.data.path = preprocessEvent(paras)
#> [1] "[2025-06-21 10:46:12.681954] Processing raw data..."
#> [1] "Input: /disk/share/lvxuan/SCSES_test//rds/"
#> [1] "Output: /disk/share/lvxuan/SCSES_test//rds_processed/"
#> [1] "Before filtering"
#> [1] "expr: 5152*15"
#> [1] "psi: 2130*15"
#> [1] "rc: 4859*15"
#> [1] "After filtering"
#> [1] "expr: 1884*15"
#> [1] "psi: 786*15"
#> [1] "rc: 1895*15"
#> [1] "[2025-06-21 10:46:13.195329] Successfully processed data."
print(rawrds.path)
#> [1] "/disk/share/lvxuan/SCSES_test//rds/"
print(processed.data.path)
#> [1] "/disk/share/lvxuan/SCSES_test//rds_processed/"

# Examine processed data
psi_processed <- readRDS(file.path(processed.data.path, "psi.rds"))
rc_processed <- readRDS(file.path(processed.data.path, "rc.rds"))

cat("Processed PSI matrix:", dim(psi_processed), "\n")
#> Processed PSI matrix: 786 15
cat("Processed RC matrix:", dim(rc_processed), "\n")
#> Processed RC matrix: 1895 15
```

Raw read count matrix, PSI matrix, and event annotation will be saved to
`work_path/rds/`.

Then, data after quality control process will be saved to
`work_path/rds_processed/`, which will be used for **subsequent
calculations**.

#### Step 5. Constructs similarity networks

To overcome the high dropout rate and limited read coverage of scRNA-seq
techniques, SCSES constructs cell similarity and event similarity
networks by K-nearest neighbor algorithm (KNN) to learn information from
similar cells and events.

``` r
# Cell similarity networks
cellnet.path = getCellSimilarity(paras)
#> [1] "[2025-06-21 10:46:13.247895] Calculate cell similarity..."
#> [1] "Input: /disk/share/lvxuan/SCSES_test//rds_processed/"
#> [1] "Output: /disk/share/lvxuan/SCSES_test//imputation/cell_similarity/"
#> [1] "feature_num=1000  checked"
#> [1] "cell_similarity_data=EXP_RBP  checked"
#> [1] "distance_method=euclidean  checked"
#> [1] "alpha_cell=0.8  checked"
#> [1] "kcell_max=8  checked"
#> [1] "kcell_min=3  checked"
#> [1] "decay_cell=0.05  checked"
#> [1] "Checking data: EXP_RBP"
#> [1] "134 rbps are used to calculate cell similarity"
#> [1] "[2025-06-21 10:46:13.301823] Computing cell similarity based on EXP_RBP"
#> [1] "Calculate similarity among 15 cells."
#> [1] "The number of features is greater than the number of rows in the input data."
#> [1] "Total 134 features will be used"
#> Delta: [0.19279533 0.03687512]
#> [1] "[2025-06-21 10:46:43.595692] Calculate cell similarity Finish."
# Event similarity networks  
eventnet.path = getEventSimilarity(paras)
#> [1] "[2025-06-21 10:46:43.602822] Calculate events KNN..."
#> [1] "Input: /disk/share/lvxuan/SCSES_test//rds_processed/"
#> [1] "Output: /disk/share/lvxuan/SCSES_test//imputation/event_similarity/"
#> [1] "alpha_event=0.8  checked"
#> [1] "kevent=5  checked"
#> [1] "decay_event=0.05  checked"
#> [1] "Checking data..."
#> [1] "Checking events..."
#> [1] "event_type=A3SS;A5SS;MXE;RI;SE  checked"
#> [1] "[2025-06-21 10:47:14.606468] Calculate events feature..."
#> [1] "[2025-06-21 10:47:14.607387] step1 Creating BSgenome for hg19======="
#> Warning: previous export ''hg19'' is being replaced
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
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> 
#> Attaching package: 'rtracklayer'
#> The following object is masked from 'package:BiocIO':
#> 
#>     FileForFormat
#> [1] "[2025-06-21 10:47:15.159689] step2 Extracting features ======="
#> [1] "[2025-06-21 10:47:15.160183] Extracting A3SS features..."
#> [1] "[2025-06-21 10:47:15.160415] Loading events..."
#> [1] "[2025-06-21 10:47:40.906637] Parsing events region..."
#> [1] "[2025-06-21 10:48:33.698315] Extracting length features."
#> [1] "[2025-06-21 10:48:33.81657] Extracting motif features."
#> [1] "[2025-06-21 10:48:33.831861] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2025-06-21 10:48:35.657475] Extracting kmer features."
#> [1] "[2025-06-21 10:48:35.700645] Extracting A Ratio features."
#> [1] "[2025-06-21 10:49:32.600886] Saving Result"
#> [1] "[2025-06-21 10:49:32.658264] Extracting A3SS features Finished"
#> [1] "[2025-06-21 10:49:32.659152] Extracting A5SS features..."
#> [1] "[2025-06-21 10:49:32.659365] Loading events..."
#> [1] "[2025-06-21 10:49:59.557972] Parsing events region..."
#> [1] "[2025-06-21 10:51:01.632601] Extracting length features."
#> [1] "[2025-06-21 10:51:01.709626] Extracting motif features."
#> [1] "[2025-06-21 10:51:01.727007] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2025-06-21 10:51:03.063393] Extracting kmer features."
#> [1] "[2025-06-21 10:51:03.078054] Extracting A Ratio features."
#> [1] "[2025-06-21 10:51:59.009568] Saving Result"
#> [1] "[2025-06-21 10:51:59.055757] Extracting A5SS features Finished"
#> [1] "[2025-06-21 10:51:59.056795] Extracting MXE features..."
#> [1] "[2025-06-21 10:51:59.057226] Loading events..."
#> [1] "[2025-06-21 10:52:27.289582] Parsing events region..."
#> [1] "[2025-06-21 10:53:25.824466] Extracting length features."
#> [1] "[2025-06-21 10:53:25.886165] Extracting motif features."
#> [1] "[2025-06-21 10:53:25.919215] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2025-06-21 10:53:27.103045] Extracting kmer features."
#> [1] "[2025-06-21 10:53:27.119756] Extracting A Ratio features."
#> [1] "[2025-06-21 10:53:27.119961] Saving Result"
#> [1] "[2025-06-21 10:53:27.125623] Extracting MXE features Finished"
#> [1] "[2025-06-21 10:53:27.125858] Extracting RI features..."
#> [1] "[2025-06-21 10:53:27.126015] Loading events..."
#> [1] "[2025-06-21 10:53:45.861901] Parsing events region..."
#> [1] "[2025-06-21 10:54:46.933083] Extracting length features."
#> [1] "[2025-06-21 10:54:47.164056] Extracting motif features."
#> [1] "[2025-06-21 10:54:47.174732] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2025-06-21 10:54:50.04962] Extracting kmer features."
#> [1] "[2025-06-21 10:54:50.075939] Extracting A Ratio features."
#> [1] "[2025-06-21 10:55:16.944149] Saving Result"
#> [1] "[2025-06-21 10:55:17.063064] Extracting RI features Finished"
#> [1] "[2025-06-21 10:55:17.063729] Extracting SE features..."
#> [1] "[2025-06-21 10:55:17.063964] Loading events..."
#> [1] "[2025-06-21 10:55:45.027923] Parsing events region..."
#> [1] "[2025-06-21 10:57:02.202212] Extracting length features."
#> [1] "[2025-06-21 10:57:03.43174] Extracting motif features."
#> [1] "[2025-06-21 10:57:03.460047] Extracting conservation features."
#> [1] "Checking chromosome prefix..."
#> [1] "[2025-06-21 10:57:13.119067] Extracting kmer features."
#> [1] "[2025-06-21 10:57:13.157136] Extracting A Ratio features."
#> [1] "[2025-06-21 10:57:13.157481] Saving Result"
#> [1] "[2025-06-21 10:57:13.430198] Extracting SE features Finished"
#> [1] "[2025-06-21 10:57:13.430578] step3 Combining events feature ======="
#> [1] "[2025-06-21 10:57:13.430785] Parsing A3SS features..."
#> [1] "[2025-06-21 10:57:13.541589] Parsing A3SS features Finished"
#> [1] "[2025-06-21 10:57:13.542378] Parsing A5SS features..."
#> [1] "[2025-06-21 10:57:13.613811] Parsing A5SS features Finished"
#> [1] "[2025-06-21 10:57:13.614898] Parsing MXE features..."
#> [1] "[2025-06-21 10:57:13.65279] Parsing MXE features Finished"
#> [1] "[2025-06-21 10:57:13.653684] Parsing RI features..."
#> [1] "[2025-06-21 10:57:13.860466] Parsing RI features Finished"
#> [1] "[2025-06-21 10:57:13.860842] Parsing SE features..."
#> [1] "[2025-06-21 10:57:14.452354] Parsing SE features Finished"
#> [1] "[2025-06-21 10:57:14.454955] step4 Encoding events feature ======="
#> [1] "[2025-06-21 10:57:14.455861] A3SS event encoding..."
#> 1/2 [==============>...............] - ETA: 0s2/2 [==============================] - 0s 2ms/step
#> [1] "[2025-06-21 10:57:25.517522] A3SS event encoding Finish."
#> [1] "[2025-06-21 10:57:25.518329] A5SS event encoding..."
#> 1/1 [==============================] - ETA: 0s1/1 [==============================] - 0s 82ms/step
#> [1] "[2025-06-21 10:57:33.311689] A5SS event encoding Finish."
#> [1] "[2025-06-21 10:57:33.312219] MXE event encoding..."
#> 1/1 [==============================] - ETA: 0s1/1 [==============================] - 0s 75ms/step
#> [1] "[2025-06-21 10:57:41.15544] MXE event encoding Finish."
#> [1] "[2025-06-21 10:57:41.156186] RI event encoding..."
#> 1/4 [======>.......................] - ETA: 0s4/4 [==============================] - 0s 2ms/step
#> [1] "[2025-06-21 10:57:52.898465] RI event encoding Finish."
#> [1] "[2025-06-21 10:57:52.899267] SE event encoding..."
#>  1/19 [>.............................] - ETA: 1s19/19 [==============================] - 0s 2ms/step
#> [1] "[2025-06-21 10:58:03.528219] SE event encoding Finish."
#> [1] "[2025-06-21 10:58:03.529057] step5 Calculate splicing regulation distance and Combine distance ======="
#> [1] "134 rbps are used to calculate splicing regulation information"
#> [1] "Save data"
#> [1] "Save data Finished"
#> [1] "[2025-06-21 10:59:26.820986] step6 Calculate combined event similarity ======="
#> [1] "[2025-06-21 10:59:29.603231] Calculate  A3SS event Similarity"
#> 
#> Attaching package: 'Matrix'
#> The following object is masked from 'package:S4Vectors':
#> 
#>     expand
#> [1] "[2025-06-21 10:59:43.621615] Calculate A3SS event Similarity Finished"
#> [1] "[2025-06-21 10:59:45.52404] Calculate  A5SS event Similarity"
#> [1] "[2025-06-21 10:59:58.664551] Calculate A5SS event Similarity Finished"
#> [1] "[2025-06-21 11:00:00.582568] Calculate  MXE event Similarity"
#> [1] "[2025-06-21 11:00:13.137459] Calculate MXE event Similarity Finished"
#> [1] "[2025-06-21 11:00:15.122637] Calculate  RI event Similarity"
#> [1] "[2025-06-21 11:00:28.183296] Calculate RI event Similarity Finished"
#> [1] "[2025-06-21 11:00:29.952301] Calculate  SE event Similarity"
#> [1] "[2025-06-21 11:00:42.973985] Calculate SE event Similarity Finished"
#> [1] "[2025-06-21 11:00:42.974376] Calculate event similarity Finished."

# Examine network structure
print(cellnet.path)
#> [1] "/disk/share/lvxuan/SCSES_test//imputation/cell_similarity/"
cell_networks <- readRDS(file.path(cellnet.path, "cell.similars.rds"))
event_networks <- readRDS(file.path(eventnet.path, "event.similars.rds"))
## three different features can be chosen to quantify cell similarity, including raw event PSI(PSI), and raw junction read counts(RC), and RBP expression(EXP_RBP)
cat("Cell similarity types:", names(cell_networks), "\n")
#> Cell similarity types: EXP_RBP
## different splicing event types
cat("Event types:", names(event_networks), "\n")
#> Event types: A3SS A5SS MXE RI SE
```

A list of event similarity for different types of splicing events will
be saved to `work_path/imputation/event_similarity/event.similars.rds`;

A list of different types of cell similarity and a list of the number of
neighbors will be saved to rds file
to`work_path/imputation/cell_similarity/cell.similars.rds` and
`work_path/imputation/cell_similarity/dyk.cell.rds`

In this step, some **parameters** can be adjusted in config file or the
function parameters directly:

For **cell similarity networks**, SCSES can use RBP expressions, Raw
read count or Raw PSI to measure cell similarities.

| Parameter | Default | Description |
|----|----|----|
| `feature_num` | 1000 | Number of high-variable features for PCA |
| `cell_similarity_data` | “EXP_RBP” | Data types for cell similarity |
| `distance_method` | “euclidean” | Distance metric |
| `alpha_cell` | 0.8 | 1 - Random walk restart probability |
| `kcell_min` | 5 | Minimum cell neighbors |
| `kcell_max` | 50 | Maximum cell neighbors |
| `decay_cell` | 0.05 | Convergence threshold |

For **event similarity networks**, Event similarities are defined by the
RBP regulatory correlations and an embedding representation by
integrating event sequence similarities.

| Parameter     | Default | Description                         |
|---------------|---------|-------------------------------------|
| `kevent`      | 10      | Number of event neighbors           |
| `alpha_event` | 0.8     | 1 - Random walk restart probability |
| `decay_event` | 0.05    | Convergence threshold               |

#### Step 6. Imputation

Based on these weighted similarity networks, SCSES next will use three
imputation strategies to aggregate the information across similar cells
or events to impute read count or PSI value

``` r
# Three-strategy imputation
Imputed.data.path = ImputationAll(paras)
#> [1] "[2025-06-21 11:00:44.883927] Get imputed result using cell similarity and event similarity."
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
#> [1] "cell_similarity_data=EXP_RBP  checked"
#> Warning in dir.create(output_path): '/disk/share/lvxuan/SCSES_test//imputation'
#> already exists
#> [1] "Output: /disk/share/lvxuan/SCSES_test//imputation/"
#> [1] "[2025-06-21 11:00:46.852427] Running Event_type=A3SS;cell_similarity_feature=EXP_RBP"
#> [1] "[2025-06-21 11:00:46.857765] Save data"
#> [1] "[2025-06-21 11:00:46.875521] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//imputation//imputation_data_3243192-500000288.h5 /disk/share/lvxuan/SCSES_test//imputation//imputation_result_3243192-500000288.mat >> /disk/share/lvxuan/SCSES_test//imputation//mat_Imputation.log 2>&1"
#> [1] "[2025-06-21 11:01:10.399819] Running Event_type=A5SS;cell_similarity_feature=EXP_RBP"
#> [1] "[2025-06-21 11:01:10.402676] Save data"
#> [1] "[2025-06-21 11:01:10.421468] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//imputation//imputation_data_3243192-500004466.h5 /disk/share/lvxuan/SCSES_test//imputation//imputation_result_3243192-500004466.mat >> /disk/share/lvxuan/SCSES_test//imputation//mat_Imputation.log 2>&1"
#> [1] "[2025-06-21 11:01:27.845692] Running Event_type=MXE;cell_similarity_feature=EXP_RBP"
#> [1] "[2025-06-21 11:01:27.846703] Save data"
#> [1] "[2025-06-21 11:01:27.860995] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//imputation//imputation_data_3243192-499996103.h5 /disk/share/lvxuan/SCSES_test//imputation//imputation_result_3243192-499996103.mat >> /disk/share/lvxuan/SCSES_test//imputation//mat_Imputation.log 2>&1"
#> [1] "[2025-06-21 11:01:49.569937] Running Event_type=RI;cell_similarity_feature=EXP_RBP"
#> [1] "[2025-06-21 11:01:49.572466] Save data"
#> [1] "[2025-06-21 11:01:49.591659] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//imputation//imputation_data_3243192-500002224.h5 /disk/share/lvxuan/SCSES_test//imputation//imputation_result_3243192-500002224.mat >> /disk/share/lvxuan/SCSES_test//imputation//mat_Imputation.log 2>&1"
#> [1] "[2025-06-21 11:02:12.385537] Running Event_type=SE;cell_similarity_feature=EXP_RBP"
#> [1] "[2025-06-21 11:02:12.393032] Save data"
#> [1] "[2025-06-21 11:02:12.420549] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/scses/run_scses.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//imputation//imputation_data_3243192-499958234.h5 /disk/share/lvxuan/SCSES_test//imputation//imputation_result_3243192-499958234.mat >> /disk/share/lvxuan/SCSES_test//imputation//mat_Imputation.log 2>&1"
#> [1] "[2025-06-21 11:02:31.272822] Get imputed result using cell similarity and event similarity Finish."

# Examine imputation results
Imputed_seperated = readRDS(Imputed.data.path)
str(Imputed_seperated,max.level=3)
#> List of 2
#>  $ cell      :List of 2
#>   ..$ EXP_RBP_PSI: num [1:786, 1:15] 0.915 0.471 0.574 0.876 0.473 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ EXP_RBP_RC : num [1:786, 1:15] 1 0.535 0.817 0.71 1 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>  $ cell_event:List of 1
#>   ..$ EXP_RBP_PSI: num [1:786, 1:15] 0.544 0.6 0.543 0.53 0.551 ...
#>   .. ..- attr(*, "dimnames")=List of 2
# explain:
# For example:
# cell/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw PSI with cell similarities.
# cell/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw inclusion and exclusion read counts with cell similarities, and then calculate the imputed PSI.
# cell_event/PSI_PSI: PSI value is used to quantify cell-cell splicing similarity. Impute raw inclusion and exclusion read counts with cell similarities, and then calculate the imputed PSI. Further impute the results using event similarities.
```

Results of each imputation strategy will be saved
into`work_path/imputation`.

#### Step 7. Estimation

##### 7.1. Estimation based on pre-trained model

We recommend different imputation strategies for four scenarios defined
by the abundance of reads counts in the target cell and neighbor cells
(ND, BD, TD+Info, and TD-Info). SCSES pre-trains models to predict the
probability of specific scenario for each cell-event pair. Finally,
SCSES calculates the PSI value using a linear combination of predictions
from the four strategies, weighted by these probabilities.

``` r
# Final estimation with scenario classification
# rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
Imputed.data.final.path = Estimation(paras,rds_imputed_file = Imputed.data.path)
#> [1] "[2025-06-21 11:02:31.352066] Combine imputed psi."
#> [1] "[2025-06-21 11:02:31.352485] Loading data..."
#> [1] "Input: /disk/share/lvxuan/SCSES_test//imputation//Imputed_seperated_499996563.rds"
#> [1] "Checking data..."
#> [1] "rc checked"
#> [1] "psi checked"
#> [1] "expr checked"
#> [1] "event checked"
#> [1] "cell similarity checked"
#> [1] "dynamic cell knn checked"
#> [1] "Pre-trained model will be used."
#> [1] "classifer checked"
#> [1] "Output: /disk/share/lvxuan/SCSES_test//imputation/"
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=EXP_RBP  checked"
#> [1] "[2025-06-21 11:02:33.508113] Combine imputed psi Finish."

# Examine final results
Imputed_combined = readRDS(Imputed.data.final.path)
str(Imputed_combined,max.level=2)
#> List of 1
#>  $ EXP_RBP: num [1:786, 1:15] 1 0.471 0.574 0.876 1 ...
#>   ..- attr(*, "dimnames")=List of 2
# The finally imputed PSI were named by cell-cell splicing similarity features, including raw event PSI(PSI), and raw junction read counts(RC), and RBP expression(EXP_RBP).
```

A list of final imputation of PSI values will be saved by serialized R
object format (.rds) in `work_path/imputation/Imputed_combined*`, where
**\*** is a random number representing different execution. The `.rds`
file can be loaded in R environment by `readRDS` function.

##### 7.2. Estimation based on fine-tuned model

To improve the fitness of models for a new dataset, we also provide a
procedure to fine-tune the model. For this analysis, we first build a
reference using a set of splicing events with conserved splicing levels
in different human tissues (<https://zenodo.org/records/6408906>). Then
we compare the splicing level in a new dataset with the reference
records, and give the scenarios definition to each event-cell pair,
which is used to fine-tune the pre-trained model.

The commands to perform these analyses:

``` r
# Final estimation with scenario classification
ftrc.path = getFtRawRC(paras)
#> [1] "Loading splicing events for classifer fine tune..."
#> [1] "Checking cells..."
#> [1] "15 cells are considered."
#> [1] "Output: /disk/share/lvxuan/SCSES_test//splicing_value_ft/"
#> [1] "[2025-06-21 11:02:33.541655] Counting reads of A3SS events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value_ft//A3SS_rjm' already exists
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
#> [1] "[2025-06-21 11:03:03.526749] Counting reads of A3SS events Finish."
#> [1] "[2025-06-21 11:03:03.527112] Counting reads of A5SS events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value_ft//A5SS_rjm' already exists
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
#> [1] "[2025-06-21 11:03:33.589083] Counting reads of A5SS events Finish."
#> [1] "[2025-06-21 11:03:33.589473] Counting reads of MXE events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value_ft//MXE_rjm' already exists
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
#> [1] "[2025-06-21 11:04:01.860919] Counting reads of MXE events Finish."
#> [1] "[2025-06-21 11:04:01.86135] Counting reads of SE events..."
#> Warning in dir.create(outpath_per_cell):
#> '/disk/share/lvxuan/SCSES_test//splicing_value_ft//SE_rjm' already exists
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
#> [1] "[2025-06-21 11:04:30.597954] Counting reads of SE events Finish."
ftpsi.path = getFtRawPSI(paras)
#> [1] "Checking raw reads..."
#> [1] "Loading splicing events for classifer fine tune..."
#> [1] "[2025-06-21 11:04:30.609568] Calculating PSI value of A3SS events..."
#> [1] "[2025-06-21 11:04:30.612852] Calculating PSI value of A3SS events Finish."
#> [1] "[2025-06-21 11:04:30.613102] Calculating PSI value of A5SS events..."
#> [1] "[2025-06-21 11:04:30.629161] Calculating PSI value of A5SS events Finish."
#> [1] "[2025-06-21 11:04:30.629678] Calculating PSI value of MXE events..."
#> [1] "[2025-06-21 11:04:30.631177] Calculating PSI value of MXE events Finish."
#> [1] "[2025-06-21 11:04:30.631451] Calculating PSI value of SE events..."
#> [1] "[2025-06-21 11:04:30.635654] Calculating PSI value of SE events Finish."
ftrds.path = mergeFtSplicingValue(paras)
ftmodel.path = FtClassifier(paras)
#> [1] "Reading true Ft PSI..."
#> [1] "Loading Pre-training classifer..."
#> [1] "[2025-06-21 11:04:30.723374] Classifer fine tune"
#> [1] "[2025-06-21 11:04:30.723963] Processing raw Ft data..."
#> [1] "Checking data..."
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=EXP_RBP  checked"
#> [1] "Calculating Classifier Features..."
#> [1] "[2025-06-21 11:04:32.993363] Save data"
#> [1] "[2025-06-21 11:04:33.009583] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//classifer//imputation1_data_3243192-500012479.h5 /disk/share/lvxuan/SCSES_test//classifer//imputation1_result_3243192-500012479.mat >> /disk/share/lvxuan/SCSES_test//classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2025-06-21 11:05:10.870352] Save data"
#> [1] "[2025-06-21 11:05:10.887563] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//classifer//imputation1_data_3243192-500021749.h5 /disk/share/lvxuan/SCSES_test//classifer//imputation1_result_3243192-500021749.mat >> /disk/share/lvxuan/SCSES_test//classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2025-06-21 11:05:48.430875] Save data"
#> [1] "[2025-06-21 11:05:48.448767] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//classifer//imputation1_data_3243192-500024016.h5 /disk/share/lvxuan/SCSES_test//classifer//imputation1_result_3243192-500024016.mat >> /disk/share/lvxuan/SCSES_test//classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2025-06-21 11:06:24.36141] Save data"
#> [1] "[2025-06-21 11:06:24.375075] Save data Finished"
#> [1] "bash /tmp/RtmpysvzUv/temp_libpath31782430f4c0d0/SCSES/matlab/imputation1/run_imputation1.sh /disk/software/matlab2022 /disk/share/lvxuan/SCSES_test//classifer//imputation1_data_3243192-499998536.h5 /disk/share/lvxuan/SCSES_test//classifer//imputation1_result_3243192-499998536.mat >> /disk/share/lvxuan/SCSES_test//classifer//mat_calculateFtFeature.log 2>&1"
#> [1] "[2025-06-21 11:07:01.843674]  Model training;similarity_type=EXP_RBP"
#> [1] "[2025-06-21 11:07:03.448754] Classifer fine tune Finish."
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
ImputedFt.data.final.path = Estimation(paras,rds_imputed_file = Imputed.data.path)
#> [1] "[2025-06-21 11:07:03.450292] Combine imputed psi."
#> [1] "[2025-06-21 11:07:03.450469] Loading data..."
#> [1] "Input: /disk/share/lvxuan/SCSES_test//imputation//Imputed_seperated_499996563.rds"
#> [1] "Checking data..."
#> [1] "rc checked"
#> [1] "psi checked"
#> [1] "expr checked"
#> [1] "event checked"
#> [1] "cell similarity checked"
#> [1] "dynamic cell knn checked"
#> [1] "Fine tune model will be used."
#> [1] "classifer checked"
#> [1] "Output: /disk/share/lvxuan/SCSES_test//imputation/"
#> [1] "Checking cell similarity type"
#> [1] "cell_similarity_data=EXP_RBP  checked"
#> [1] "[2025-06-21 11:07:05.461887] Combine imputed psi Finish."

print(ImputedFt.data.final.path)
#> [1] "/disk/share/lvxuan/SCSES_test//imputation//Imputed_combined_500023390.rds"

# Examine final results
Imputed_combined = readRDS(ImputedFt.data.final.path)
str(Imputed_combined,max.level=2)
#> List of 1
#>  $ EXP_RBP: num [1:786, 1:15] 1 0.471 0.574 0.876 1 ...
#>   ..- attr(*, "dimnames")=List of 2
Imputed_combined[["EXP_RBP"]][1:3,1:3]
#>                                                                                                                                                                                                                                         SRR11826368.bam
#> isoform1=exon:chr1:151139802-151139890:+@junction:chr1:151139891-151140624:+@exon:chr1:151140625-151140814:+|isoform2=exon:chr1:151139802-151139890:+@junction:chr1:151139891-151140619:+@exon:chr1:151140620-151140814:+|SCNM1|A3SS          1.0000000
#> isoform1=exon:chr1:153610771-153610924:+@junction:chr1:153610925-153614721:+@exon:chr1:153614722-153614905:+|isoform2=exon:chr1:153610771-153610924:+@junction:chr1:153610925-153614718:+@exon:chr1:153614719-153614905:+|CHTOP|A3SS          0.4709326
#> isoform1=exon:chr1:169772310-169772450:+@junction:chr1:169772451-169773252:+@exon:chr1:169773253-169773381:+|isoform2=exon:chr1:169772310-169772450:+@junction:chr1:169772451-169773215:+@exon:chr1:169773216-169773381:+|C1orf112|A3SS       0.5737308
#>                                                                                                                                                                                                                                         SRR11826371.bam
#> isoform1=exon:chr1:151139802-151139890:+@junction:chr1:151139891-151140624:+@exon:chr1:151140625-151140814:+|isoform2=exon:chr1:151139802-151139890:+@junction:chr1:151139891-151140619:+@exon:chr1:151140620-151140814:+|SCNM1|A3SS          1.0000000
#> isoform1=exon:chr1:153610771-153610924:+@junction:chr1:153610925-153614721:+@exon:chr1:153614722-153614905:+|isoform2=exon:chr1:153610771-153610924:+@junction:chr1:153610925-153614718:+@exon:chr1:153614719-153614905:+|CHTOP|A3SS          0.4473823
#> isoform1=exon:chr1:169772310-169772450:+@junction:chr1:169772451-169773252:+@exon:chr1:169773253-169773381:+|isoform2=exon:chr1:169772310-169772450:+@junction:chr1:169772451-169773215:+@exon:chr1:169773216-169773381:+|C1orf112|A3SS       1.0000000
#>                                                                                                                                                                                                                                         SRR11826409.bam
#> isoform1=exon:chr1:151139802-151139890:+@junction:chr1:151139891-151140624:+@exon:chr1:151140625-151140814:+|isoform2=exon:chr1:151139802-151139890:+@junction:chr1:151139891-151140619:+@exon:chr1:151140620-151140814:+|SCNM1|A3SS          1.0000000
#> isoform1=exon:chr1:153610771-153610924:+@junction:chr1:153610925-153614721:+@exon:chr1:153614722-153614905:+|isoform2=exon:chr1:153610771-153610924:+@junction:chr1:153610925-153614718:+@exon:chr1:153614719-153614905:+|CHTOP|A3SS          0.4941710
#> isoform1=exon:chr1:169772310-169772450:+@junction:chr1:169772451-169773252:+@exon:chr1:169773253-169773381:+|isoform2=exon:chr1:169772310-169772450:+@junction:chr1:169772451-169773215:+@exon:chr1:169773216-169773381:+|C1orf112|A3SS       0.8639692
```

#### Step 8. Cell Clustering

Here is an example of UMAP visualization based on test data:

``` r
install.packages("umap")
```

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

# Load data and annotations

if(!file.exists(paste0(system.file("analysis", package = "SCSES"),'/annotation.txt'))) {
  stop("Annotation file not found")
}

if(!file.exists(ImputedFt.data.final.path)) {
  stop("Imputed data file not found")
}

annotation=read.table(paste0(system.file("analysis", package = "SCSES"),'/annotation.txt'),sep="\t")
# Imputed_combined=readRDS(paste0(system.file("analysis", package = "SCSES"),'/Imputed_combined_499997048.rds'))
Imputed_combined=readRDS(ImputedFt.data.final.path)

data_umap=calcu_umap(Imputed_combined[[1]],n_neighbors = 3)
row.names(data_umap)=gsub(".bam","",row.names(data_umap))
data_umap$group=annotation$V2[match(row.names(data_umap),annotation$V1)]

p=ggplot(data = data_umap,aes(x =V1 ,y =V2))+
  geom_point(aes(fill=group),shape=21,size=1.5,stroke=0.05)+
  scale_fill_manual(values = mycol)+
  xlab("UMAP1")+
  ylab("UMAP2")

print(p)
```

<img src="man/figures/README-Cell Clustering-1.png" width="100%" />
