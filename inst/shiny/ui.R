library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(parallel)
source("www/R/uiFunction.R")

header <- dashboardHeader(title = "SCSES")
sidebar <- dashboardSidebar(
       sidebarMenu(
              menuItem(text = "Configuration", selected = T, tabName = "configuration")
       )
)
body <- dashboardBody(
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
                                          textInput(inputId = "dataset", label = "Dataset", placeholder = "Input the dataset name"),
                                          fluidRow(
                                                 createTextFile(id = "bam_path", label = "Bam File Path", type = "shinyDirButton"),
                                                 column(
                                                        width = 2,
                                                        numericInput(inputId = "readlength", label = "Read Length", value = 125, min = 0)
                                                 ),
                                                 column(
                                                        width = 2,
                                                        tags$label(class = "control-label", "placehold", "for" = "read.paired", style = "visibility: hidden")
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "read.paired", label = "<b>Paired Read</b>", value = T, onLabel = "paired", offLabel = "single", labelWidth = "150px")
                                                 ),
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "read.umi", label = "<b>Full-length</b>", value = T, onLabel = "full_length", offLabel = "UMI", labelWidth = "100px")
                                                 )
                                          )
                                   ),
                                   tabPanel(
                                          "Tools Info",
                                          h3("SCSES Program Essential"),
                                          textInput(inputId = "conda_envname", label = "conda environment", placeholder = "Input the name of conda env for SCSES",value='base'),
                                          fluidRow(
                                                 createTextFile(id = "condabin_path", label = "Conda bin Path", type = "shinyDirButton",value="/software/miniconda3/bin/"),
                                                 createTextFile(id = "python_path", label = "Python Path", type = "shinyFilesButton"),
                                                 createTextFile(id = "JAVA_path", label = "JAVA Path", type = "shinyFilesButton"),
                                                 createTextFile(id = "MCR_path", label = "MCR Path", type = "shinyDirButton",value='/opt/mcr/R2022b/')
                                          ),
                                          fluidRow(
                                                 createTextFile(id = "Samtools_path", label = "Samtools Path", type = "shinyFilesButton"),
                                                 createTextFile(id = "featurecounts_path", label = "FeatureCounts Path", type = "shinyFilesButton")
                                          ),
                                          hr(),
                                          h3("Splicing Event Detection Essential"),
                                          fluidRow(
                                                 createTextFile(id = "rmats_path", label = "rMats Path", type = "shinyFilesButton",value='/software/rmats_turbo_v4_3_0/rmats.py')
                                          ),
                                          textInput(inputId = "MAJIQ_env", label = "MAJIQ conda environment", placeholder = "Input the name of conda env for MAJIQ",value = 'MAJIQ'),
                                          fluidRow(
                                                 createTextFile(id = "MAJIQ_license_path", label = "MAJIQ License Path", type = "shinyFilesButton",value='/software/majiq_license_academic_official.lic')
                                          ),
                                          fluidRow(
                                                 createTextFile(id = "IRFinder_path", label = "IRFinder Path", type = "shinyFilesButton"),
                                                 createTextFile(id = "STAR_path", label = "STAR Path", type = "shinyFilesButton")
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
                                                 createTextFile(id = "work_path", label = "Work Path", type = "shinyDirButton"),
                                                 column(
                                                        width = 2,
                                                        numericInput(inputId = "max.core", label = "Maximal CPU Cores", value = 1, min = 1, max = detectCores())
                                                 )
                                          ),
                                          hr(),
                                          h3("Events Filter in Merged Bam"),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "exon.intron.read", label = "Minimal Read Count on Exon-Intron Boundary", value = 150, min = 0)
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "junction.read", label = "Minimal Read Count on Junction", value = 150, min = 0)
                                                 )
                                          ),
                                          hr(),
                                          h3("Events Filter in Single Bam"),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "exp.dropout.ratio", label = "Maximal Dropout Ratio of Genes Expression", value = 0.9, min = 0, max = 1, step = 0.01)
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "psi.dropout.ratio", label = "Maximal Dropout Ratio of Splicing Events", value = 0.9, min = 0, max = 1, step = 0.01)
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "cell.read", label = "Minimal Read Count", value = 1000, min = 0)
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "cell.gene.exp", label = "Minimal Count of Expressed Gene", value = 1, min = 0)
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "cell.mt.pct", label = "Maximal Ratio of Expressed MT-Genes", value = 0.5, min = 0, max = 1, step = 0.01)
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "minRC", label = "Minimal Read Count Supporting an Event", value = 5, min = 0)
                                                 ),
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "minCell", label = "Minimal Cell Count Containing an Event", value = 20, min = 0)
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "filter.mt", label = "<b>Filter Out Mitochondrial Genes</b>", value = T, onLabel = "TRUE", offLabel = "FALSE", labelWidth = "150px")
                                                 ),
                                                 column(
                                                        width = 3,
                                                        switchInput(inputId = "filter.rp", label = "<b>Filter Out ribosomal Genes</b>", value = T, onLabel = "TRUE", offLabel = "FALSE", labelWidth = "150px")
                                                 )
                                          )
                                   ),
                                   tabPanel(
                                          "Imputation",
                                          h3("Reference Files"),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        textInput(inputId = "ref_name", label = "Reference Name", value = "hg19")
                                                 ),
                                                 createTextFile(id = "fa_path", label = "Reference Genome(.fa)", type = "shinyFilesButton"),
                                                 createTextFile(id = "GTF_path", label = "GTF File", type = "shinyFilesButton"),
                                                 createTextFile(id = "GFF_path", label = "GFF File", type = "shinyFilesButton"),
                                                 createTextFile(id = "RBP_path", label = "RBP File", type = "shinyFilesButton")
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
                                                        )
                                                 )
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 4,
                                                        numericInput(inputId = "feature_num", label = "The Number of High Variable Cell Similarity Features", value = 1000, min = 10)
                                                 )
                                          ),
                                          hr(),
                                          h3("Event Similarity"),
                                          fluidRow(
                                                 column(
                                                        width = 12,
                                                        checkboxGroupButtons(
                                                               inputId = "event.types", label = "Event Types",
                                                               choices = c("SE" = "SE", "RI" = "RI", "A3SS" = "A3SS", "A5SS" = "A5SS", "MXE" = "MXE", "AL" = "AL"), selected = c("SE", "RI", "A3SS", "A5SS", "MXE"),
                                                               checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))
                                                        )
                                                 )
                                          ),
                                          fluidRow(
                                                 createTextFile(id = "phast_path", label = "PhastCons Path", type = "shinyFilesButton")
                                          ),
                                          fluidRow(
                                                 column(
                                                        width = 5,
                                                        switchInput(inputId = "remove_chr", label = "<b>Remove chromosome prefix</b>", value = T, onLabel = "TRUE", offLabel = "FALSE", labelWidth = "180px")
                                                 ),
                                                 column(
                                                   width = 5,
                                                   switchInput(inputId = "chr_prefix", label = "<b>Add chr prefix when extracting conservation scores</b>", value = T, onLabel = "TRUE", offLabel = "FALSE", labelWidth = "180px")
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
                                                                      numericInput(inputId = "SE.epoch", label = "Epoch", value = 100, min = 50),
                                                                      numericInput(inputId = "SE.embedding", label = "Embedding Dim", value = 32, min = 50),
                                                                      textInput(inputId = "SE.layer", label = "Layers", value = "256,128")
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
                                                                      numericInput(inputId = "RI.epoch", label = "Epoch", value = 100, min = 50),
                                                                      numericInput(inputId = "RI.embedding", label = "Embedding Dim", value = 32, min = 50),
                                                                      textInput(inputId = "RI.layer", label = "Layers", value = "256,128")
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
                                                                      numericInput(inputId = "A3SS.epoch", label = "Epoch", value = 100, min = 50),
                                                                      numericInput(inputId = "A3SS.embedding", label = "Embedding Dim", value = 32, min = 50),
                                                                      textInput(inputId = "A3SS.layer", label = "Layers", value = "256,128")
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
                                                                      numericInput(inputId = "A5SS.epoch", label = "Epoch", value = 100, min = 50),
                                                                      numericInput(inputId = "A5SS.embedding", label = "Embedding Dim", value = 32, min = 50),
                                                                      textInput(inputId = "A5SS.layer", label = "Layers", value = "256,128")
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
                                                       numericInput(inputId = "MXE.epoch", label = "Epoch", value = 100, min = 50),
                                                       numericInput(inputId = "MXE.embedding", label = "Embedding Dim", value = 32, min = 50),
                                                       textInput(inputId = "MXE.layer", label = "Layers", value = "256,128")
                                                     )
                                                   )
                                                 ),
                                                 column(
                                                        width = 2,
                                                        tags$fieldset(
                                                               style = "border: 1px solid;border-color:#C0C0C0",
                                                               tags$legend("AL", style = "border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';", align = "center"),
                                                               div(
                                                                      style = "margin-left:1rem;margin-right:1rem",
                                                                      numericInput(inputId = "AL.epoch", label = "Epoch", value = 100, min = 50),
                                                                      numericInput(inputId = "AL.embedding", label = "Embedding Dim", value = 32, min = 50),
                                                                      textInput(inputId = "AL.layer", label = "Layers", value = "256,128")
                                                               )
                                                        )
                                                 )
                                          ),
                                          hr(),
                                          h3("Network Fusion Parameters"),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "all.decay", label = "Decay", value = 0.05, min = 0)
                                                 )
                                          ),
                                          h4("Cell Similarity Network"),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        radioGroupButtons(
                                                               inputId = "similar.method", label = "Similarity Type",
                                                               choices = c("Euclidean" = "euclidean", "Cosine" = "cosine"), selected = "euclidean", status = "default", individual = T,
                                                               checkIcon = list(
                                                                      yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                                                                      no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                                                               )
                                                        )
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "cell.kmin", label = "Min K cell", value = 5, min = 0, step = 1)
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "cell.kmax", label = "Max K cell", value = 20, min = 0, step = 1)
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "cell.alpha", label = "Cell Alpha", value = 0.8, min = 0, max = 1, step = 0.1)
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "cell.decay", label = "Cell Decay", value = 0.05, min = 0)
                                                 )
                                          ),
                                          h4("Event Similarity Network"),
                                          fluidRow(
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "event.k", label = "K Neighbor Event", value = 5, min = 0, step = 1)
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "event.alpha", label = "Cell Alpha", value = 0.8, min = 0, max = 1, step = 0.1)
                                                 ),
                                                 column(
                                                        width = 3,
                                                        numericInput(inputId = "event.decay", label = "Event Decay", value = 0.05, min = 0)
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
