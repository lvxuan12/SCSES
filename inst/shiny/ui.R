library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(parallel)
source('www/R/uiFunction.R')

header=dashboardHeader(title = 'SCESC')
sidebar=dashboardSidebar(
  sidebarMenu(
    menuItem(text = 'Configuration',selected = T,tabName = 'configuration'),
    menuItem(text = 'Excutaton',selected = F,tabName = 'excute')
  )
)
body=dashboardBody(
  tabItems(
    tabItem(tabName = 'configuration',
            fluidPage(tags$head(tags$script(src = 'js/shinyEvents.js')),
              tabBox(id = 'basic.info',title = 'Basic Configuration',width = 12,
                     tabPanel("Data Info",
                              h3('Data Information'),
                              textInput(inputId='dataset',label='Dataset',placeholder = 'Input the dataset name'),
                              fluidRow(
                                       createTextFile(id = 'cellinfo_file',label = 'Cell Info File'),
                                       div(class='col-sm-9',
                                           radioGroupButtons(inputId = 'cellinfo.column',label = 'Cell Name Column',choices = NA,
                                                             checkIcon = list(yes = icon("ok",lib = "glyphicon"),no = icon("remove",lib = "glyphicon"))
                                          )
                                      )
                              ),
                              fluidRow(
                                createTextFile(id = 'bam_path',label = 'Bam File Path',type="shinyDirButton"),
                                column(width=2,
                                       numericInput(inputId = 'readlength',label = 'Read Length', value = 125,min=0)
                                ),
                                column(width=2,
                                       tags$label(class="control-label","placehold","for"='read.paired',style="visibility: hidden")
                                )
                              ),
                              fluidRow(
                                column(width=3,
                                      switchInput(inputId = 'read.paired',label = '<b>Paired Read</b>',value = T,onLabel = 'Paired',offLabel = 'Single',labelWidth = '150px')
                                ),
                                column(width=2,
                                        switchInput(inputId = 'read.umi',label = '<b>Full/UMI</b>',value = T,onLabel = 'Full',offLabel = 'UMI',labelWidth = '100px')
                                )
                              )
                              
                     ),
                     tabPanel("Tools Info",
                              h3('SCESC Program Essential'),
                              fluidRow(
                                createTextFile(id = 'SCESC_src',label = 'SCESC Source Directory',type="shinyDirButton"),
                                createTextFile(id = 'Rscript_path',label = 'Rscript Path',type="shinyFilesButton"),
                                createTextFile(id = 'JAVA_path',label = 'JAVA Path',type="shinyFilesButton"),
                                createTextFile(id = 'Samtools_path',label = 'Samtools Path',type="shinyFilesButton"),
                                createTextFile(id = 'featurecounts_path',label = 'FeatureCounts Path',type="shinyFilesButton")
                              ),
                              fluidRow(
                                # createTextFile(id = 'Matlab_path',label = 'Matlab Path',type="shinyFilesButton"),
                                createTextFile(id = 'MCR_path',label = 'MCR Path',type="shinyDirButton")
                              ),
                              hr(),
                              h3("Splicing Event Detection Essential"),
                              fluidRow(
                                createTextFile(id = 'rmats_path',label = 'rMats Path',type="shinyFilesButton"),
                                createTextFile(id = 'MAJIQ_path',label = 'MAJIQ Path',type="shinyFilesButton"),
                                createTextFile(id = 'IRFinder_path',label = 'IRFinder Path',type="shinyFilesButton")
                              )
                     ),
                     tabPanel("Others",
                              fluidRow(
                                createTextFile(id = 'BRIE_path',label = 'BRIE Path',type="shinyFilesButton"),
                                createTextFile(id = 'outrigger_path',label = 'Outrigger Path',type="shinyFilesButton")
                              )
                     )
              )
            ),
            fluidPage(
              tabBox(id = 'task.info',title = 'Task Configuration',width = 12,
                     tabPanel("Basic",
                              h3('Reference Files'),
                              fluidRow(
                                createTextFile(id = 'work_path',label = 'Work Path',type="shinyDirButton"),
                                column(width=2,
                                       numericInput(inputId='max.core',label = 'Maximal CPU Cores',value = 1,min = 1,max = detectCores())  
                                )
                              ),
                              hr(),
                              h3('Events Filter in Merged Bam'),
                              fluidRow(
                                column(width=3,
                                       numericInput(inputId='exon.intron.read',label = 'Minimal Read Count on Exon-Intron Boundary',value = 150,min = 0)
                                ),
                                column(width=3,
                                       numericInput(inputId='junction.read',label = 'Minimal Read Count on Junction',value = 150,min = 0)
                                       )
                              ),
                              hr(),
                              h3('Events Filter in Single Bam'),
                              fluidRow(
                                column(width=3,
                                       numericInput(inputId='exp.dropout.ratio',label = 'Maximal Dropout Ratio of Genes Expression',value = 0.99,min = 0,max = 1,step = 0.01)
                                ),
                                column(width=3,
                                       numericInput(inputId='psi.dropout.ratio',label = 'Maximal Dropout Ratio of Splicing Events',value = 0.99,min = 0,max = 1,step = 0.01)
                                )
                              ),
                              fluidRow(
                                column(width=3,
                                       numericInput(inputId='cell.read',label = 'Minimal Read Count per Cell',value = 1000,min = 0)
                                ),
                                column(width=3,
                                       numericInput(inputId='cell.gene.exp',label = 'Minimal Count of Expressed Gene per Cell',value = 1,min = 0)
                                ),
                                column(width=3,
                                       numericInput(inputId='cell.mt.pct',label = 'Maximal Ratio of Expressed MT-Genes per Cell',value = 0.5,min = 0,max = 1,step = 0.01)
                                )
                              )
                     ),
                     tabPanel("Imputation",
                              h3('Reference Files'),
                              fluidRow(
                                column(width=3,
                                       textInput(inputId='ref_name',label = 'Reference Name',value = 'hg19')       
                                ),
                                createTextFile(id = 'fa_path',label = 'Reference Genome(.fa)',type="shinyFilesButton"),
                                createTextFile(id = 'GTF_path',label = 'GTF File',type="shinyFilesButton"),
                                createTextFile(id = 'GFF_path',label = 'GFF File',type="shinyFilesButton"),
                                createTextFile(id = 'RBP_path',label = 'RBP File',type="shinyFilesButton"),
                                createTextFile(id = 'finetune_raw_path',label = 'Fine-tune Raw File',type="shinyFilesButton")
                              ),
                              hr(),
                              h3('Cell Similarity'),
                              fluidRow(
                                column(width=12,
                                       checkboxGroupButtons(inputId = 'cell.similar',label = 'Cell Similarity Feature',
                                                            choices = c("Gene Expression"="EXP",'PSI'='PSI','Read Count'='RC',"RBP Expression"='EXP_RBP'),selected = 'EXP_RBP',
                                                            checkIcon = list(yes = icon("ok",lib = "glyphicon"),no = icon("remove",lib = "glyphicon"))
                                       )  
                                )
                                
                              ),
                              hr(),
                              h3('Event Similarity'),
                              fluidRow(
                                column(width=12,
                                       checkboxGroupButtons(inputId = 'event.types',label = 'Event Types',
                                                     choices = c("SE"="SE",'RI'='RI','A3SS'='A3SS',"A5SS"='A5SS','MXE'='MXE','ALE'='ALE'),selected = c('SE','RI','A3SS','A5SS','MXE'),
                                                     checkIcon = list(yes = icon("ok",lib = "glyphicon"),no = icon("remove",lib = "glyphicon"))
                                       )
                                )
                              ),
                              fluidRow(
                                column(width=3,
                                       textInput(inputId='ref.pkg',label = 'Reference Package',value = 'hg19')
                                ),
                                createTextFile(id = 'phast_path',label = 'PhastCons Path',type="shinyFilesButton"),
                                column(width=3,
                                       textInput(inputId='chr.prefix',label = 'Chromosome Prefix',value = 'chr')
                                ),
                                column(width=3,
                                       div(class='form-group shiny-input-container',
                                           tags$label(class="control-label","placehold","for"="remove_chr",style="visibility: hidden"),
                                           switchInput(inputId = 'remove_chr',label = '<b>Remove Chr?</b>',value = F,onLabel = 'Yes',offLabel = 'No',labelWidth = '100px')
                                       )
                                  )
                              ),
                              h4('Autoencoder Parameters'),
                              fluidRow(
                                column(width=2,
                                       tags$fieldset(style="border: 1px solid;border-color:#C0C0C0",
                                         tags$legend("SE",style="border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';",align='center'),
                                         div(style='margin-left:1rem;margin-right:1rem',
                                           numericInput(inputId='SE.epoch',label = 'Epoch',value = 100,min = 50),
                                           numericInput(inputId='SE.embedding',label = 'Embedding Dim',value = 100,min = 50),
                                           textInput(inputId='SE.layer',label = 'Layers',value = '256,128')  
                                         ) 
                                       )
                                ),
                                column(width=2,
                                       tags$fieldset(style="border: 1px solid;border-color:#C0C0C0",
                                                     tags$legend("RI",style="border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';",align='center'),
                                                     div(style='margin-left:1rem;margin-right:1rem',
                                                         numericInput(inputId='RI.epoch',label = 'Epoch',value = 100,min = 50),
                                                         numericInput(inputId='RI.embedding',label = 'Embedding Dim',value = 100,min = 50),
                                                         textInput(inputId='RI.layer',label = 'Layers',value = '256,128')  
                                                     ) 
                                       )
                                ),
                                column(width=2,
                                       tags$fieldset(style="border: 1px solid;border-color:#C0C0C0",
                                                     tags$legend("A3SS",style="border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';",align='center'),
                                                     div(style='margin-left:1rem;margin-right:1rem',
                                                         numericInput(inputId='A3SS.epoch',label = 'Epoch',value = 100,min = 50),
                                                         numericInput(inputId='A3SS.embedding',label = 'Embedding Dim',value = 100,min = 50),
                                                         textInput(inputId='A3SS.layer',label = 'Layers',value = '256,128')  
                                                     ) 
                                       )
                                ),
                                column(width=2,
                                       tags$fieldset(style="border: 1px solid;border-color:#C0C0C0",
                                                     tags$legend("A5SS",style="border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';",align='center'),
                                                     div(style='margin-left:1rem;margin-right:1rem',
                                                         numericInput(inputId='A5SS.epoch',label = 'Epoch',value = 100,min = 50),
                                                         numericInput(inputId='A5SS.embedding',label = 'Embedding Dim',value = 100,min = 50),
                                                         textInput(inputId='A5SS.layer',label = 'Layers',value = '256,128')  
                                                     ) 
                                       )
                                ),
                                column(width=2,
                                       tags$fieldset(style="border: 1px solid;border-color:#C0C0C0",
                                                     tags$legend("MXE",style="border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';",align='center'),
                                                     div(style='margin-left:1rem;margin-right:1rem',
                                                         numericInput(inputId='MXE.epoch',label = 'Epoch',value = 100,min = 50),
                                                         numericInput(inputId='MEX.embedding',label = 'Embedding Dim',value = 100,min = 50),
                                                         textInput(inputId='MEX.layer',label = 'Layers',value = '256,128')  
                                                     ) 
                                       )
                                ),
                                column(width=2,
                                       tags$fieldset(style="border: 1px solid;border-color:#C0C0C0",
                                          tags$legend("ALE",style="border: 1px #C0C0C0 solid;background-color:#fff;color:#000;width:50px;text-align='center';",align='center'),
                                          div(style='margin-left:1rem;margin-right:1rem',
                                             numericInput(inputId='AL.epoch',label = 'Epoch',value = 100,min = 50),
                                             numericInput(inputId='AL.embedding',label = 'Embedding Dim',value = 100,min = 50),
                                             textInput(inputId='AL.layer',label = 'Layers',value = '256,128')  
                                         )
                                      )
                              )
                        ),
                        hr(),
                        h3('Network Fusion Parameters'),
                        fluidRow(
                          column(width=3,
                                  radioGroupButtons(inputId = 'dim.reduction',label = 'Dimension Reduction Method',
                                                    choices = c('PCA'='PCA','tSNE'='tSNE','UMAP'='UMAP'),selected = 'PCA',status = 'default',individual = T,
                                                    checkIcon = list(yes = tags$i(class = "fa fa-circle",style = "color: steelblue"),
                                                                     no = tags$i(class = "fa fa-circle-o",style = "color: steelblue")))
                          ),
                          column(width=3,
                                 radioGroupButtons(inputId = 'similar.method',label = 'Similarity Type',
                                                   choices = c('Euclidean'='euclidean'),selected = 'euclidean',status = 'default',individual = T,
                                                   checkIcon = list(yes = tags$i(class = "fa fa-circle",style = "color: steelblue"),
                                                                    no = tags$i(class = "fa fa-circle-o",style = "color: steelblue")))
                          ),
                          column(width=3,
                                 numericInput(inputId='all.decay',label = 'Decay',value = 0.05,min = 0)
                          )
                        ),
                        h4('Cell Similarity Network'),
                        fluidRow(
                          column(width=3,
                                 numericInput(inputId='cell.kmin',label = 'Min K cell',value = 5,min = 0,step=1)
                          ),
                          column(width=3,
                                 numericInput(inputId='cell.kmax',label = 'Max K cell',value = 5,min = 0,step=1)    
                          ),
                          column(width=3,
                               numericInput(inputId='cell.alpha',label = 'Cell Alpha',value = 0.8,min = 0,max = 1,step = 0.1)
                          ),
                          column(width=3,
                                 numericInput(inputId='cell.decay',label = 'Cell Decay',value = 0.05,min = 0)
                          )
                        ),
                        h4('Event Similarity Network'),
                        fluidRow(
                          column(width=3,
                                 numericInput(inputId='event.k',label = 'K Neighbor Event',value = 10,min = 0,step=1)
                          ),
                          column(width=3,
                                 numericInput(inputId='event.alpha',label = 'Cell Alpha',value = 0.8,min = 0,max = 1,step = 0.1)
                          ),
                          column(width=3,
                                 numericInput(inputId='event.decay',label = 'Event Decay',value = 0.05,min = 0)
                          )
                        )
                    ),
                     tabPanel(title='Others',
                              h3('BRIE Setting'),
                              fluidRow(
                                createTextFile(id = 'brie_ref',label = 'BRIE Reference',type = 'shinyFilesButton'),
                                createTextFile(id = 'brie_feature',label = 'BRIE Feature',type = 'shinyFilesButton')
                              ),
                              hr(),
                              h3('Outrigger Setting'),
                              fluidRow(
                                createTextFile(id = 'outrigger_size',label = 'Genome Size',type = 'shinyFilesButton')
                              )
                              
                     )
              )
            ),
            fluidPage(
              fluidRow(
                column(width=4,
                       actionBttn(inputId = 'create_configure',label = 'Create Config',color = 'primary')    
                )
                
              )
            )
            
    ),
    tabItem(tabName = 'excute')
  )
)

ui <- dashboardPage(
  header = header,
  sidebar = sidebar,
  body = body,
  title = 'SCESC'
)
