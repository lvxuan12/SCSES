library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(parallel)
library(prompter)
library(dplyr)
source("www/R/uiFunction.R")

header <- dashboardHeader(title = "SCSES")
sidebar <- dashboardSidebar(
       sidebarMenu(
              menuItem(text = "Configuration", selected = T, tabName = "configuration")
       )
)
body <- dashboardBody(
       use_prompt(),
       tabItems(
              tabItem(
                     tabName = "configuration",
                     fluidPage(
                            tags$head(tags$script(src = "js/shinyEvents.js")),
                            tabBox(
                                   id = "basic.info", title = "Basic Configuration", width = 12,
                                   tabPanel(
                                          "Data Info",
                                          h3("Data Information"),
                                          textInput(inputId = "dataset", label = "Dataset", placeholder = "Input the dataset name")%>%
                                            add_prompt(position = 'bottom',message = 'Dataset name, used for naming the configuration file and the gene expression file output by featureCount.',type = 'info',size='large'),

                                          fluidRow(
                                                 createTextFile(id = "bam_path", label = "Directory containing Bam files",msg='Directory to bam files. refer to "Basic/bam_path" in configure file.', type = "shinyDirButton"),
                                                 column(
                                                        width = 2,
                                                        numericInput(inputId = "readlength", label = "Read Length", value = 100, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'The sequence length of cells. refer to "Baisc/readlength" in configure file.',type = 'info',size='large')

                                                 ),
                                                 column(
                                                        width = 2,
                                                        tags$label(class = "control-label", "placehold", "for" = "read.paired", style = "visibility: hidden")
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "read.paired", label = "<b>Paired Read</b>", value = T, onLabel = "paired", offLabel = "single", labelWidth = "150px")%>%
                                                          add_prompt(position = 'bottom',message = 'Sequence type: "paired" for paired-end, "single" for single-end. Refer to "Baisc/readlength" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "read.umi", label = "<b>Full-length</b>", value = T, onLabel = "full_length", offLabel = "UMI", labelWidth = "100px")%>%
                                                          add_prompt(position = 'bottom',message = 'Sequence technique: full-length (e.g. SMART-seq) or UMI-based (e.g. 10x genomics). Refer to "Baisc/sequence" in configure file.',type = 'info',size='large')
                                                 )
                                          )
                                   ),
                                   tabPanel(
                                          "Tools Info",
                                          h3("SCSES Program Essential"),
                                          textInput(inputId = "conda_envname", label = "conda environment", value = 'base',placeholder = "Input the name of conda env for SCSES")%>%
                                            add_prompt(position = 'bottom',message = 'The conda environment used for running SCSES,default is "base". Refer to "Baisc/conda_envname" in configure file.',type = 'info',size='large'),
                                          fluidRow(
                                                 createTextFile(id = "condabin_path", label = "Conda bin Path", type = "shinyDirButton",msg = 'Full path of "conda/bin" directory. Refer to "Basic/conda_binpath" in configure file.'),
                                                 createTextFile(id = "python_path", label = "Python Path", type = "shinyFilesButton",msg='Full path of "python" program. It will be filled by the result of "which python" automatically. Refer to "Basic/python_path" in configure file.'),
                                                 createTextFile(id = "JAVA_path", label = "JAVA Path", type = "shinyFilesButton",msg='Full path of "JAVA" program. It will be filled by the result of "which java" automatically. Refer to "Basic/java_path" in configure file.'),
                                                 createTextFile(id = "MCR_path", label = "MCR Path", type = "shinyDirButton",msg = 'Full path of Matlab compiler runtime. Refer to "Basic/mcr_path" in the configure file.')
                                          ),
                                          fluidRow(
                                                 createTextFile(id = "Samtools_path", label = "Samtools Path", type = "shinyFilesButton",msg='Full path of "samtools" program. It will be filled by the result of "which samtools" automatically. Refer to "Basic/samtools" in configure file.'),
                                                 createTextFile(id = "featurecounts_path", label = "FeatureCounts Path", type = "shinyFilesButton",'Full path of "featureCounts" program. It will be filled by the result of "which featureCounts" automatically. Refer to "Basic/featureCounts_path" in configure file.')
                                          ),
                                          hr(),
                                          h3("Splicing Event Detection Essential"),
                                          fluidRow(
                                                 createTextFile(id = "rmats_path", label = "rMats Path", type = "shinyFilesButton",msg = 'Full path of "rMATS" program. It will be filled by the result of "which rmats" automatically. Refer to "Basic/rMATS_path" in configure file.')
                                          ),
                                          textInput(inputId = "MAJIQ_env", label = "MAJIQ conda environment", placeholder = "Input the name of conda env for MAJIQ")%>%
                                            add_prompt(position = 'bottom',message = 'The conda environment name for running MAJIQ. Refer to "Basic/MAJIQ_env" in configure file.',type = 'info',size='large'),
                                          fluidRow(
                                                 createTextFile(id = "MAJIQ_license_path", label = "MAJIQ License Path", type = "shinyFilesButton",msg='The full path of MAJIQ license file. Refer to "Basic/majiq_license_file" in configure file.')
                                          ),
                                          fluidRow(
                                                 createTextFile(id = "IRFinder_path", label = "IRFinder Path", type = "shinyFilesButton",msg='Full path of "IRFinder" program. It will be filled by the result of "which IRFinder" automatically. Refer to "Basic/IRFinder_path" in configure file.'),
                                                 createTextFile(id = "STAR_path", label = "STAR Path", type = "shinyFilesButton",msg='Full path of "STAR" program. It will be filled by the result of "which STAR" automatically. Refer to "Basic/STAR_path" in configure file.')
                                          )
                                   )
                            )
                     ),
                     fluidPage(
                            tabBox(
                                   id = "task.info", title = "Task Configuration", width = 12,
                                   tabPanel(
                                          "Basic",
                                          fluidRow(
                                                 createTextFile(id = "work_path", label = "Work Path", type = "shinyDirButton",msg='Working directory, where all the software outputs will be stored. Refer to "Basic/work_path" in configure file.'),
                                                 column(
                                                        width = 2,
                                                        numericInput(inputId = "max.core", label = "Maximal CPU Cores", value = 1, min = 1, max = detectCores())%>%
                                                          add_prompt(position = 'bottom',message = 'Thread count for parallel processing. Refer to "Basic/core" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          hr(),
                                          h3("Events Filter in Merged Bam"),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "exon.intron.read", label = "Minimal Read Count on Exon-Intron Boundary", value = 25, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'The minimum number of total reads spanning exon and intron regions required to a RI event. Refer to "Basic/filter_merged_bam/ExonToIntronReads" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "junction.read", label = "Minimal Read Count on Junction", value = 25, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'The minimum number of total junction reads required to support a non-RI event. Refer to "Basic/filter_merged_bam/junctionReads" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          hr(),
                                          h3("Events Filter in Single Bam"),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "exp.dropout.ratio", label = "Maximal Dropout Ratio of Genes Expression", value = 0.9, min = 0, max = 1, step = 0.01)%>%
                                                          add_prompt(position = 'bottom',message = 'Maximum allowable percentage of gene dropouts. Refer to "Basic/filter_sc/min.percentCells.gene" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "psi.dropout.ratio", label = "Maximal Dropout Ratio of Splicing Events", value = 0.9, min = 0, max = 1, step = 0.01)%>%
                                                          add_prompt(position = 'bottom',message = 'Maximum allowable percentage of event dropouts (PSI=0,1,or NA). Refer to "Basic/filter_sc/min.percentCells.event" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "cell.read", label = "Minimal Read Count", value = 1, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'Minimum read counts per junction required in at least minCell cells. Refer to "Basic/filter_sc/min.nCount" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "cell.gene.exp", label = "Minimal Count of Expressed Gene", value = 1, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'Minimum read count required to consider a gene expressed in a cell. Refer to "Basic/filter_sc/min.nFeatures" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "cell.mt.pct", label = "Maximal Ratio of Expressed MT-Genes", value = 1, min = 0, max = 1, step = 0.01)%>%
                                                          add_prompt(position = 'bottom',message = 'Maximum allowable fraction of mitochondrial gene counts per cell. Refer to "Basic/filter_sc/max.percentMT" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "minRC", label = "Minimal Read Count Supporting an Event", value = 5, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'Minimum read counts per junction required in at least minCell cells. Refer to "Basic/filter_sc/max.percentMT" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "minCell", label = "Minimal Cell Count Containing an Event", value = 5, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'Minimum number of cells per junction required to have at least min.RC reads. Refer to "Basic/filter_sc/minCell" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "filter.mt", label = "<b>Remove MT Genes</b>", value = T, onLabel = "TRUE", offLabel = "FALSE", labelWidth = "150px")%>%
                                                          add_prompt(position = 'bottom',message = 'If remove mitochondrial genes or not. Refer to "Basic/filter_sc/filter.mt" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "filter.rp", label = "<b>Remove Ribo Genes</b>", value = T, onLabel = "TRUE", offLabel = "FALSE", labelWidth = "150px")%>%
                                                          add_prompt(position = 'bottom',message = 'If remove ribosomal genes or not. Refer to "Basic/filter_sc/filter.rp" in configure file.',type = 'info',size='large')
                                                 )
                                          )
                                   ),
                                   tabPanel(
                                          "Imputation",
                                          h3("Reference Files"),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        textInput(inputId = "ref_name", label = "Reference Name", value = "hg19")%>%
                                                          add_prompt(position = 'bottom',message = 'Genome name (e.g. hg19/hg38/mm10) required in following situations:\n
                                                                     (1) running MAJIQ; \n
                                                                     (2) naming a forge BSgenome package for sequence feature extraction; \n
                                                                     (3) deciding the event coordition used for model fine-tuning (only support hg19/hg38/mm10 in this part).
                                                                     Refer to "Basic/refgenome/genome_name" in configure file.',type = 'info',size='large')
                                                 ),
                                                 createTextFile(id = "fa_path", label = "Reference Genome(.fa)", type = "shinyFilesButton",msg = 'Full path to the fasta file of reference genome. Refer to "Basic/refgenome/ref_path" in configure file.'),
                                                 createTextFile(id = "GTF_path", label = "GTF File", type = "shinyFilesButton",msg = 'Full path to the gene annotation in GTF format. Refer to "Basic/refgenome/gtf_path" in configure file.'),
                                                 createTextFile(id = "GFF_path", label = "GFF File", type = "shinyFilesButton",msg = 'Full path to the gene annotation in GFF format. Refer to "Basic/refgenome/gff_path" in configure file.'),
                                                 createTextFile(id = "RBP_path", label = "RBP File", type = "shinyFilesButton",msg = 'Full path to RBP file used for cell similarity and RBP regulatory features. Refer to "Task/impute/rbp" in configure file.')
                                          ),
                                          hr(),
                                          h3("Cell Similarity"),
                                          fluidRow(
                                                 column(
                                                        width = 12,
                                                        checkboxGroupButtons(
                                                               inputId = "cell.similar", label = "Cell Similarity Feature",
                                                               choices = c("PSI" = "PSI", "Read Count" = "RC", "RBP Expression" = "EXP_RBP"), selected = "EXP_RBP",
                                                               checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))
                                                        )%>%
                                                          add_prompt(position = 'bottom',message = 'The feature type used to calculate cell similarity. Refer to "Task/impute/cell_similarity_data" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "feature_num", label = "Highly Variable Feature Count", value = 1000, min = 10)%>%
                                                          add_prompt(position = 'bottom',message = 'The number of highly variable features used for PCA. Refer to "Task/impute/feature_num" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          hr(),
                                          h3("Event Similarity"),
                                          fluidRow(
                                                 column(
                                                        width = 12,
                                                        checkboxGroupButtons(
                                                               inputId = "event.types", label = "Event Types",
                                                               choices = c("SE" = "SE", "RI" = "RI", "A3SS" = "A3SS", "A5SS" = "A5SS", "MXE" = "MXE"), selected = c("SE", "RI", "A3SS", "A5SS", "MXE"),
                                                               checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))
                                                        )%>%
                                                          add_prompt(position = 'bottom',message = 'Splicing event types for detection and imputation. Refer to "Task/event/event_type" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          fluidRow(
                                                 createTextFile(id = "phast_path", label = "PhastCons Path", type = "shinyFilesButton",msg = 'Full path of sequence conservertion file (phastCons.bw). Refer to "Task/impute/event_features/phast_path" in configure file')
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "remove_chr", label = '<b>Remove "chr" prefix</b>', value = F, onLabel = "FALSE", offLabel = "TRUE", labelWidth = "180px")%>%
                                                          add_prompt(position = 'bottom',message = 'If removing "chr" prefix from chromosome numbers. Refer to "Task/event/event_type" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                   width = 3,
                                                   switchInput(inputId = "chr_prefix", label = '<b>Add "chr" prefix</b>', value = F, onLabel = "FALSE", offLabel = "TRUE", labelWidth = "180px")%>%
                                                     add_prompt(position = 'bottom',message = 'If add "chr" as prefix to event id to ensure the chromosome id in bam files, splicing events and sequence conservation file to be the same. Refer to "Task/impute/event_features/chr_prefix" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          #textInput(inputId = "chr_prefix", label = "Add chromosome prefix", placeholder = paste("Add chr prefix when extracting conservation scores")),
                                          h4("Autoencoder Parameters"),
                                          fluidRow(
                                                 column(
                                                        width = 2,
                                                        tags$fieldset(
                                                               style = "border: 1px solid;border-color:#C0C0C0",
                                                               tags$legend("SE", style = "border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';", align = "center"),
                                                               div(
                                                                      style = "margin-left:1rem;margin-right:1rem",
                                                                      numericInput(inputId = "SE.epoch", label = "Epoch", value = 100, min = 50)%>%
                                                                        add_prompt(position = 'bottom',message = 'Training epoch for SE events. Refer to "Task/impute/event_features/AE/SE/epoch" in configure file.',type = 'info',size='large'),
                                                                      numericInput(inputId = "SE.embedding", label = "Embedding Dim", value = 32, min = 50)%>%
                                                                        add_prompt(position = 'bottom',message = 'Latent embedding dimention for SE events. Refer to "Task/impute/event_features/AE/SE/embedding" in configure file.',type = 'info',size='large'),
                                                                      textInput(inputId = "SE.layer", label = "Layers", value = "256,128")%>%
                                                                        add_prompt(position = 'bottom',message = 'Dimention of each layer in Encoder for SE events, format:layer1,layer2,...,layerN. Refer to "Task/impute/event_features/AE/SE/layer" in configure file.',type = 'info',size='large')
                                                               )
                                                        )
                                                 ),
                                                 column(
                                                        width = 2,
                                                        tags$fieldset(
                                                               style = "border: 1px solid;border-color:#C0C0C0",
                                                               tags$legend("RI", style = "border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';", align = "center"),
                                                               div(
                                                                      style = "margin-left:1rem;margin-right:1rem",
                                                                      numericInput(inputId = "RI.epoch", label = "Epoch", value = 100, min = 50)%>%
                                                                        add_prompt(position = 'bottom',message = 'Training epoch for RI events. Refer to "Task/impute/event_features/AE/RI/epoch" in configure file.',type = 'info',size='large'),
                                                                      numericInput(inputId = "RI.embedding", label = "Embedding Dim", value = 32, min = 50)%>%
                                                                        add_prompt(position = 'bottom',message = 'Latent embedding dimention for RI events. Refer to "Task/impute/event_features/AE/RI/embedding" in configure file.',type = 'info',size='large'),
                                                                      textInput(inputId = "RI.layer", label = "Layers", value = "256,128")%>%
                                                                        add_prompt(position = 'bottom',message = 'Dimention of each layer in Encoder for RI events, format:layer1,layer2,...,layerN. Refer to "Task/impute/event_features/AE/RI/layer" in configure file.',type = 'info',size='large')
                                                               )
                                                        )
                                                 ),
                                                 column(
                                                        width = 2,
                                                        tags$fieldset(
                                                               style = "border: 1px solid;border-color:#C0C0C0",
                                                               tags$legend("A3SS", style = "border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';", align = "center"),
                                                               div(
                                                                      style = "margin-left:1rem;margin-right:1rem",
                                                                      numericInput(inputId = "A3SS.epoch", label = "Epoch", value = 100, min = 50)%>%
                                                                        add_prompt(position = 'bottom',message = 'Training epoch for A3SS events. Refer to "Task/impute/event_features/AE/A3SS/epoch" in configure file.',type = 'info',size='large'),
                                                                      numericInput(inputId = "A3SS.embedding", label = "Embedding Dim", value = 32, min = 50)%>%
                                                                        add_prompt(position = 'bottom',message = 'Latent embedding dimention for A3SS events. Refer to "Task/impute/event_features/AE/A3SS/embedding" in configure file.',type = 'info',size='large'),
                                                                      textInput(inputId = "A3SS.layer", label = "Layers", value = "256,128")%>%
                                                                        add_prompt(position = 'bottom',message = 'Dimention of each layer in Encoder for A3SS events, format:layer1,layer2,...,layerN. Refer to "Task/impute/event_features/AE/A3SS/layer" in configure file.',type = 'info',size='large')
                                                               )
                                                        )
                                                 ),
                                                 column(
                                                        width = 2,
                                                        tags$fieldset(
                                                               style = "border: 1px solid;border-color:#C0C0C0",
                                                               tags$legend("A5SS", style = "border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';", align = "center"),
                                                               div(
                                                                      style = "margin-left:1rem;margin-right:1rem",
                                                                      numericInput(inputId = "A5SS.epoch", label = "Epoch", value = 100, min = 50)%>%
                                                                        add_prompt(position = 'bottom',message = 'Training epoch for A5SS events. Refer to "Task/impute/event_features/AE/A5SS/epoch" in configure file.',type = 'info',size='large'),
                                                                      numericInput(inputId = "A5SS.embedding", label = "Embedding Dim", value = 32, min = 50)%>%
                                                                        add_prompt(position = 'bottom',message = 'Latent embedding dimention for A5SS events. Refer to "Task/impute/event_features/AE/A5SS/embedding" in configure file.',type = 'info',size='large'),
                                                                      textInput(inputId = "A5SS.layer", label = "Layers", value = "256,128")%>%
                                                                        add_prompt(position = 'bottom',message = 'Dimention of each layer in Encoder for A5SS events, format:layer1,layer2,...,layerN. Refer to "Task/impute/event_features/AE/A5SS/layer" in configure file.',type = 'info',size='large')
                                                               )
                                                        )
                                                 ),
                                                 column(
                                                   width = 2,
                                                   tags$fieldset(
                                                     style = "border: 1px solid;border-color:#C0C0C0",
                                                     tags$legend("MXE", style = "border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';", align = "center"),
                                                     div(
                                                       style = "margin-left:1rem;margin-right:1rem",
                                                       numericInput(inputId = "MXE.epoch", label = "Epoch", value = 100, min = 50)%>%
                                                         add_prompt(position = 'bottom',message = 'Training epoch for MXE events. Refer to "Task/impute/event_features/AE/MXE/epoch" in configure file.',type = 'info',size='large'),
                                                       numericInput(inputId = "MXE.embedding", label = "Embedding Dim", value = 32, min = 50)%>%
                                                         add_prompt(position = 'bottom',message = 'Latent embedding dimention for MXE events. Refer to "Task/impute/event_features/AE/MXE/embedding" in configure file.',type = 'info',size='large'),
                                                       textInput(inputId = "MXE.layer", label = "Layers", value = "256,128")%>%
                                                         add_prompt(position = 'bottom',message = 'Dimention of each layer in Encoder for MXE events, format:layer1,layer2,...,layerN. Refer to "Task/impute/event_features/AE/MXE/layer" in configure file.',type = 'info',size='large')
                                                     )
                                                   )
                                                 )#,
                                                 # column(
                                                 #        width = 2,
                                                 #        tags$fieldset(
                                                 #               style = "border: 1px solid;border-color:#C0C0C0",
                                                 #               tags$legend("AL", style = "border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';", align = "center"),
                                                 #               div(
                                                 #                      style = "margin-left:1rem;margin-right:1rem",
                                                 #                      numericInput(inputId = "AL.epoch", label = "Epoch", value = 100, min = 50),
                                                 #                      numericInput(inputId = "AL.embedding", label = "Embedding Dim", value = 32, min = 50),
                                                 #                      textInput(inputId = "AL.layer", label = "Layers", value = "256,128")
                                                 #               )
                                                 #        )
                                                 # )
                                          ),
                                          hr(),
                                          h3("Network Fusion Parameters"),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "all.decay", label = "Decay", value = 0.05, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'Convergence threshold of imputation. Refer to "Task/impute/decay_impute" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          h4("Cell Similarity Network"),
                                          fluidRow(
                                                 # column(
                                                 #        width = 3,
                                                 #        radioGroupButtons(
                                                 #               inputId = "similar.method", label = "Similarity Type",
                                                 #               choices = c("Euclidean" = "euclidean", "Cosine" = "cosine"), selected = "euclidean", status = "default", individual = T,
                                                 #               checkIcon = list(
                                                 #                      yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                                                 #                      no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                                                 #               )
                                                 #        )
                                                 # ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "cell.kmin", label = "Min K cell", value = 3, min = 0, step = 1)%>%
                                                          add_prompt(position = 'bottom',message = 'Minimum number of dynamic cell neighbors. Refer to "Task/impute/KNN/cell/kmin" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "cell.kmax", label = "Max K cell", value = 8, min = 0, step = 1)%>%
                                                          add_prompt(position = 'bottom',message = 'Maximal number of dynamic cell neighbors. Refer to "Task/impute/KNN/cell/kmax" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "cell.alpha", label = "Random walk probability (Non restart)", value = 0.8, min = 0, max = 1, step = 0.1)%>%
                                                          add_prompt(position = 'bottom',message = 'Random walk probability (1-restart probability) for cell similarity network. Refer to "Task/impute/KNN/cell/alpha" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "cell.decay", label = "Cell Decay", value = 0.05, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'Convergence threshold of cell similarity network diffusion. Refer to "Task/impute/KNN/cell/decay" in configure file.',type = 'info',size='large')
                                                 )
                                          ),
                                          h4("Event Similarity Network"),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "event.k", label = "K Neighbor Event", value = 5, min = 0, step = 1)%>%
                                                          add_prompt(position = 'bottom',message = 'K value for event neighbors. Refer to "Task/impute/KNN/event/k" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "event.alpha", label = "Random walk probability (Non restart)", value = 0.8, min = 0, max = 1, step = 0.1)%>%
                                                          add_prompt(position = 'bottom',message = 'Random walk probability (1-restart probability) for event similarity network. Refer to "Task/impute/KNN/event/alpha" in configure file.',type = 'info',size='large')
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "event.decay", label = "Event Decay", value = 0.05, min = 0)%>%
                                                          add_prompt(position = 'bottom',message = 'Convergence threshold of event similarity network diffusion. Refer to "Task/impute/KNN/event/decay" in configure file.',type = 'info',size='large')
                                                 )
                                          )
                                   )
                            )
                     ),
                     fluidPage(
                            fluidRow(
                                   column(
                                          width = 4,
                                          actionBttn(inputId = "create_configure", label = "Create Config", color = "primary")
                                   )
                            )
                     )
              )
       )
)

ui <- dashboardPage(
       header = header,
       sidebar = sidebar,
       body = body,
       title = "SCSES"
)
