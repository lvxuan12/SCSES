Menu
================

- [SCSES](#scses)
  - [Installation](#installation)
    - [Dependencies and requirements](#dependencies-and-requirements)
      - [1. python module](#1-python-module)
      - [2. Softwares](#2-softwares)
    - [Installation from GitHub](#installation-from-github)
      - [Tips for some Installation
        error](#tips-for-some-installation-error)
  - [SCSES input](#scses-input)
  - [Getting started](#getting-started)
    - [Run SCSES using one command:](#run-scses-using-one-command)
    - [Run SCSES step by step (smart-seq2
      data)](#run-scses-step-by-step-smart-seq2-data)
      - [Step1. Read configure file](#step1-read-configure-file)
      - [Step2. Get gene expression](#step2-get-gene-expression)
      - [Step3. Detect splicing events](#step3-detect-splicing-events)
      - [Step4. Quantify splicing
        events](#step4-quantify-splicing-events)
      - [Step5. Constructs similarity
        networks](#step5-constructs-similarity-networks)
      - [Step6. Imputation](#step6-imputation)
      - [Step7. Estimation](#step7-estimation)
  - [Tutorials](#tutorials)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCSES

<!-- badges: start -->

<!-- badges: end -->

Single-cell Splicing Estimation based on Network Diffusion

## Installation

### Dependencies and requirements

We recommend a new conda environment to install SCSES:

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
python). MAJIQ should be built in a new environment.

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

MAJIQ will not function without providing the license file.

##### 2.3 [IRFinder](https://github.com/dgaolab/IRFinder)

``` bash
wget https://github.com/dgaolab/IRFinder/archive/refs/tags/v1.3.1.tar.gz
tar -zxvf v1.3.0.tar.gz
export PATH=/path/to/IRFinder-2.0.1/bin/:$PATH
```

[STAR](https://github.com/alexdobin/STAR) will be needed to build
IRFinder reference

##### 2.4 [samtools](https://github.com/samtools/samtools)

### Installation from GitHub

SCSES is installed directly from github using the following command:

``` r
remotes::install_github("lvxuan12/SCSES")
```

#### Tips for some Installation error

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

SCSES requires four inputs from the user

##### 1. bam files sorted by coordinate

Test bam files can be downloaded from …

##### 2. genome FASTA file and annotation GTF file

##### 3. configure file

You can use `createConfigshiny` to generate configure file for your
dataset.

``` r
createConfigshiny(host, port) 
```

A json file will be generated in the work_path you provided.

##### 4. phast conservation file in bigWig format

For human and mouse, you could download it directly from UCSC browser:
[mm10.60way.phastCons.bw](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/),
[hg38.phastCons100way.bw](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/)
and
[hg19.100way.phastCons.bw](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/).

## Getting started

### Run SCSES using one command:

### Run SCSES step by step (smart-seq2 data)

#### Step1. Read configure file

``` r
library(SCSES)
#paras_file: path to configure file generated in the previous step
paras = readSCSESconfig(paras_file)
```

#### Step2. Get gene expression

##### TPM matrix (for smart-seq2 dataset)

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

##### Normalized UMI count matrix (for UMI dataset)

You can use `get10XEXPmatrix` to generate Normalized UMI count matrix
from 10X CellRanger hdf5 file, which will save normalized UMI count to
`work_path/rds/`.

``` r
rds.path = get10XEXPmatrix(paras,expr_path,sample_name)
```

#### Step3. Detect splicing events

To define a global set of all splicing events, SCSES firstly merges all
bam files from every single cell to construct a pseudo-bulk bam file,
and identifies all types of splicing events by conventional algorithms.

###### for smart-seq2 dataset

``` r
pseudobulk.path = createPseudobulk(paras)
event.path = detectEvents(paras)
```

Different types of splicing events will be saved to `work_path/events/`,
separately.

###### for UMI dataset

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

#### Step4. Quantify splicing events

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

Raw read count matrix and PSI matrix for different types of splicing
events will be saved to `work_path/splicing_value/`, separately.

#### Step5. Constructs similarity networks

To overcome the high dropout rate and limited read coverage of scRNA-seq
techniques, SCSES constructs cell similarity and event similarity
networks by K-nearest neighbor algorithm (KNN) to learn information from
similar cells/events.

``` r
cellnet.path = getCellSimilarity(paras)
eventnet.path = getEventSimilarity(paras)
```

A list of event similarity for different types of splicing events will
be saved to `work_path/imputation/event_similarity/event.similars.rds`;

A list of different types of cell similarity and a list of the number of
neighbors will be saved to rds file
to`work_path/imputation/cell_similarity/cell.similars.rds` and
`work_path/imputation/cell_similarity/dyk.cell.rds`

In this step, some parameters can be adjusted in configure file or the
function parameters directly:

##### For cell similarity networks:

SCSES can use RBP expressions, Raw read count or Raw PSI to measure cell
similarities.

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

##### For event similarity networks:

Event similarities are defined by the RBP regulatory correlations and an
embedding representation by integrating event sequence similarities.

`ae.para`: parameters of encoding sequence features

`rbp`: expression of those RBPs will be used to calculate RBP regulatory
correlations.

`kevent`: the number of neighbors

`alpha_event`: restart probability for random walk

`decay_event`: threshold of change in the similarity matrix

#### Step6. Imputation

Based on these weighted similarity networks, SCSES next aggregates the
information across similar cells or events to impute read count or PSI
value.

``` r
Imputed.data.path = ImputationAll(paras)
```

A list of three imputation strategies result will be saved to
`work_path/imputation/`.

#### Step7. Estimation

We recommend different imputation strategies for four scenarios defined
by the abundance of reads counts in the target cell and neighbor cells
(ND, BD, TD+Info, and TD-Info). SCSES pre-trains models to predict the
probability of specific scenario for each cell-event pair. Finally,
SCSES calculates the PSI value using a linear combination of predictions
from the four strategies, weighted by these probabilities.

``` r
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
Imputed.data.path = Estimation(paras,rds_imputed_file)
```

##### Fine-tune the model

To improve the fitness of models for the new dataset, we also provide a
procedure to fine-tune the model. We collect a set of splicing events
with conserved splicing levels in different human tissue. For a new
dataset, we compare the splicing level in new data with the reference
records, and give the scenarios definition to each event-cell pair,
which is used to fine-tune the pre-trained model.

``` r
ftrc.path = getFtRawRC(paras)
ftpsi.path = getFtRawPSI(paras)
ftrds.path = mergeFtSplicingValue(paras)
ftmodel.path = FtClassifier(paras)
#rds_imputed_file: path to the list of three imputation strategies results generated in the previous step
Imputed.data.path = Estimation(paras,rds_imputed_file)
```

A list of final imputation of PSI values will be saved to
`work_path/imputation/`.

## Tutorials

For nPSC example, see the Rmarkdown tutorials in `analysis/nPSC.Rmd`
