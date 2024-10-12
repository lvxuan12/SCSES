library(shiny)
library(shinydashboard)
library(shinyFiles)
library(fs)
library(jsonlite)

server=function(input,output,session){
  source('www/R/serverFunction.R',local = T)
  
  # file.field=c('cellinfo_file','bam_path','SCESC_src','Rscript_path','JAVA_path','Samtools_path','featurecounts_path',
  #              'Matlab_path','MCR_path','rmats_path','MAJIQ_path','IRFinder_path','BRIE_path','outrigger_path','work_path',
  #              'fa_path','GTF_path','GFF_path','RBP_path','finetune_raw_path','phast_path','brie_ref','brie_feature','outrigger_size')
  volumns=c(workdir='.',Home=path_home(),getVolumes())
  env=environment()
  
  raw.values=isolate(reactiveValuesToList(input))
  
  for(ff in file.attr$id)
  {
    if(file.attr[ff,'type']=='shinyFilesButton')
    {
      shinyFileChoose(input = input, id = ff, roots=volumns,session = session,defaultRoot = 'workdir')  
    }
    else
    {
      shinyDirChoose(input = input, id = ff, roots=volumns,session = session,defaultRoot = 'workdir')  
    } 
  }
  
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'SCESC_src'),value = dirname(normalizePath('./')))
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'Rscript_path'),value = getDefaultPath(cmd = 'Rscript'))
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'JAVA_path'),value = getDefaultPath(cmd = 'java'))
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'Samtools_path'),value = getDefaultPath(cmd = 'samtools'))
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'featurecounts_path'),value = getDefaultPath(cmd = 'featureCounts'))
  # updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'Matlab_path'),value = getDefaultPath(cmd = 'Matlab_path'))
  
  # textinput action----
  observeEvent(c(input$cellinfo_file,input$cellinfofile),{
      if(showSelect('cellinfo_file'))
      {
        isolate({
          cellinfofile=input$cellinfofile
        })
        updateValue("cellinfo_file")
        if(!updateValue("cellinfofile"))
        {
          if(cellinfofile!="")
          {
            if(grepl(pattern = "\\.csv$",x = cellinfofile))
            {
              sep=','
            }
            else
            {
              sep="\t"
            }
            tryCatch({
              printMSG(paste('Read New File',cellinfofile))
              data=read.table(cellinfofile,header = T,sep = sep,check.names = F)
              columns=colnames(data)
              names(columns)=columns
              updateRadioGroupButtons(session = session,inputId = 'cellinfo.column',choices = columns,selected = columns[1])  
            },error=function(e){
              print(e)
              updateRadioGroupButtons(session = session,inputId = 'cellinfo.column',choices = 'NA',selected = 'NA')  
            })
            
          }  
        }
      }
  })
  observeEvent(c(input$bam_path,input$bam_path),{
    showSelect('bam_path')
    updateValue("bam_path")
    updateValue("bampath")
    
  })
  observeEvent(c(input$SCESC_src,input$SCESCsrc),{
    showSelect('SCESC_src')
    updateValue("SCESC_src")
    updateValue("SCESCsrc")
  })
  observeEvent(c(input$Rscript_path,input$Rscriptpath),{
    showSelect('Rscript_path')
    updateValue("Rscript_path")
    updateValue("Rscriptpath")
  })
  observeEvent(c(input$JAVA_path,input$JAVApath),{
    showSelect('JAVA_path')
    updateValue("JAVA_path")
    updateValue("JAVApath")
  })
  observeEvent(c(input$Samtools_path,input$Samtoolspath),{
    showSelect('Samtools_path')
    updateValue("Samtools_path")
    updateValue("Samtoolspath")
  })
  observeEvent(c(input$featurecounts_path,input$featurecountspath),{
    showSelect('featurecounts_path')
    updateValue("featurecounts_path")
    updateValue("featurecountspath")
  })
  observeEvent(c(input$Matlab_path,input$Matlabpath),{
    showSelect('Matlab_path')
    updateValue("Matlab_path")
    updateValue("Matlabpath")
  })
  observeEvent(c(input$MCR_path,input$MCRpath),{
    showSelect('MCR_path')
    updateValue("MCR_path")
    updateValue("MCRpath")
  })
  observeEvent(c(input$rmats_path,input$rmatspath),{
    showSelect('rmats_path')
    updateValue("rmats_path")
    updateValue("rmatspath")
  })
  observeEvent(c(input$MAJIQ_path,input$MAJIQpath),{
    showSelect('MAJIQ_path')
    updateValue("MAJIQ_path")
    updateValue("MAJIQpath")
  })
  observeEvent(c(input$IRFinder_path,input$IRFinderpath),{
    showSelect('IRFinder_path')
    updateValue("IRFinder_path")
    updateValue("IRFinderpath")
  })
  observeEvent(c(input$BRIE_path,input$BRIEpath),{
    showSelect('BRIE_path')
    updateValue("BRIE_path")
    updateValue("BRIEpath")
  })
  observeEvent(c(input$outrigger_path,input$outriggerpath),{
    showSelect('outrigger_path')
    updateValue("outrigger_path")
    updateValue("outriggerpath")
  })
  observeEvent(c(input$work_path,input$workpath),{
    showSelect('work_path')
    updateValue("work_path")
    updateValue("workpath")
  })
  observeEvent(c(input$fa_path,input$fapath),{
    showSelect('fa_path')
    updateValue("fa_path")
    updateValue("fapath")
  })
  observeEvent(c(input$GTF_path,input$GTFpath),{
    showSelect('GTF_path')
    updateValue("GTF_path")
    updateValue("GTFpath")
  })
  observeEvent(c(input$GFF_path,input$GFFpath),{
    showSelect('GFF_path')
    updateValue("GFF_path")
    updateValue("GFFpath")
  })
  observeEvent(c(input$RBP_path,input$RBPpath),{
    showSelect('RBP_path')
    updateValue("RBP_path")
    updateValue("RBPpath")
  })
  observeEvent(c(input$finetune_raw_path,input$finetunerawpath),{
    showSelect('finetune_raw_path')
    updateValue("finetune_raw_path")
    updateValue("finetunerawpath")
  })
  observeEvent(c(input$phast_path,input$phastpath),{
    showSelect('phast_path')
    updateValue("phast_path")
    updateValue("phastpath")
  })
  observeEvent(c(input$brie_ref,input$brieref),{
    showSelect('brie_ref')
    updateValue("brie_ref")
    updateValue("brieref")
  })
  observeEvent(c(input$brie_feature,input$briefeature),{
    showSelect('brie_feature')
    updateValue("brie_feature")
    updateValue("briefeature")
  })
  observeEvent(c(input$outrigger_size,input$outriggersize),{
    showSelect('outrigger_size')
    updateValue("outrigger_size")
    updateValue("outriggersize")
  })
  
  # create configuration action----
  observeEvent(input$create_configure,{
    config=list()
    flag=validation(input,session)
    if(!flag)
    {
      return()
    }
    config[['DataSet']]=input$dataset
    config[['Basic']][['cell_path']]=normalizePath(input$cellinfofile)
    config[['Basic']][['cell_column']]=input$cellinfo.column
    config[['Basic']][['bam_path']]=normalizePath(input$bampath)
    config[['Basic']][['work_path']]=normalizePath(input$workpath)
    config[['Basic']][['core']]=input$max.core
    #config[['Basic']][['memory']]=input$
    config[['Basic']][['scrpath']]=normalizePath(input$SCESCsrc)
    config[['Basic']][['script']]=normalizePath(paste0(config$Basic$scrpath,'/impute/'))
    config[['Basic']][['mcr_path']]=normalizePath(input$MCRpath)
    config[['Basic']][['matlab_script']]=normalizePath(paste0(config$Basic$scrpath,'/impute/src_matlab/'))
    config[['Basic']][['jar_path']]=normalizePath(paste0(config$Basic$scrpath,'/JAVA/ScPsiImputation'))
    config[['Basic']][['Rscript_path']]=normalizePath(input$Rscriptpath)
    config[['Basic']][['samtools_path']]=normalizePath(input$Samtoolspath)
    config[['Basic']][['MAJIQ_path']]=normalizePath(input$MAJIQpath)
    config[['Basic']][['featureCounts_path']]=normalizePath(input$featurecountspath)
    config[['Basic']][['outrigger']]=normalizePath(input$outriggerpath)
    config[['Basic']][['BRIE']]=normalizePath(input$BRIEpath)
    config[['Basic']][['rMATs']]=normalizePath(input$rmatspath)
    config[['Basic']][['IRFinder_path']]=normalizePath(input$IRFinderpath)
    config[['Basic']][['filter_merged_bam']]=list(ExonToIntronReads=input$exon.intron.read,junctionReads=input$junction.read)
    config[['Basic']][['filter_events_sc']]=list(expr_0_percent=input$exp.dropout.ratio,
                                                 psi_0_or_1_percent=input$psi.dropout.ratio,
                                                 expr_sum=input$cell.read,
                                                 gene_sum=input$cell.gene.exp,
                                                 MT_gene_pct=input$cell.mt.pct)
    config[['Task']][['impute']][['gtf_path']]=normalizePath(input$GTFpath)
    config[['Task']][['impute']][['script_data_preprocess']]=normalizePath(paste0(config$Basic$script,'/func_data_preprocess.r'))
    config[['Task']][['impute']][['script_impute_func']]=normalizePath(paste0(config$Basic$script,'/Dynamic_Kcell.py'))
    config[['Task']][['impute']][['output_path']]=normalizePath(paste0(config$Basic$work_path,'/imputation/'))
    config[['Task']][['impute']][['script_impute']]=normalizePath(paste0(config$Basic$script,'/impute_final.r'))
    config[['Task']][['impute']][['script_events_feature']]=normalizePath(paste0(config$Basic$script,'/basic.event.feature.final.r'))
    config[['Task']][['impute']][['psi_gtex_path']]=normalizePath(paste0(config$Basic$scrpath,'/data/gtex/hg19/psi/psi_gtex_select.txt'))
    config[['Task']][['impute']][['cell_similarity']]=paste(input$cell.similar,collapse = ';')
    config[['Task']][['impute']][['rbp']]=normalizePath(input$RBPpath)
    config[['Task']][['impute']][['event_features']]=list(pkg=input$ref.pkg,
                                                          phast_path=normalizePath(input$phastpath),
                                                          chr_prefix=input$chr.prefix)
    
    for(type in input$event.types)
    {
      config[['Task']][['impute']][['AE_para']][[type]]=list(epoch=input[[paste0(type,'.epoch')]],
                                                 embedding=input[[paste0(type,'.embedding')]],
                                                 layer=paste0('[',input[[paste0(type,'.layer')]],']')
                                                )
    }
    config[['Task']][['impute']][['parameter']]=list(D_reduct=input$dim.reduction,
                                                     similar_method=input$similar.method,
                                                     decay_final=input$all.decay
                                                )
    config[['Task']][['impute']][['parameter']][['cell']]=list(kmax=input$cell.kman,
                                                               kmin=input$cell.kmin,
                                                               alpha=input$cell.alpha,
                                                               decay=input$cell.decay)
    config[['Task']][['impute']][['parameter']][['event']]=list(kmax=input$event.k,
                                                               alpha=input$event.alpha,
                                                               decay=input$event.decay)
    config[['Task']][['impute']][['model1']]=normalizePath(paste0(config$Basic$script,'/model_change_nonchange.rdata'))
    config[['Task']][['impute']][['model2']]=normalizePath(paste0(config$Basic$script,'/model_change_01_change_other.rdata'))
    browser()
    config[['Task']][['MAJIQ']]=list(gff_path=normalizePath(input$GFFpath),
                                     genome_name=input$ref_name,
                                     readlength=input$readlength,
                                     script=normalizePath(paste0(input$scrpath,'/events/run_majq.sh')))
    config[['Task']][['rMats']]=list(gtf_path=normalizePath(input$GTFpath),
                                     paired=input$read.paired,
                                     readlength=input$readlength,
                                     script=normalizePath(paste0(input$scrpath,'/events/run_rmats.sh')))
    config[['Task']][['IRFinder']]=list(gtf_path=normalizePath(input$GTFpath),
                                        ref_path=normalizePath(input$fapath),
                                        readlength=input$readlength,
                                        script=normalizePath(paste0(input$scrpath,'/events/run_irfinder.sh')))
    config[['Task']][['BRIE']]=list(brie_ref=normalizePath(input$brieref),
                                    brie_feature=normalizePath(input$briefeature),
                                    script=normalizePath(paste0(input$scrpath,'/events/run_brie.sh')))
    config[['Task']][['outrigger']]=list(gtf_path=normalizePath(input$GTFpath),
                                         ref_path=normalizePath(input$fapath),
                                         genome_size=normalizePath(input$outriggersize),
                                         script=normalizePath(paste0(input$scrpath,'/events/run_outrigger.sh')))
    
    config[['Task']][['featureCounts']]=list(gtf_path=normalizePath(input$GTFpath),
                                         ref_path=normalizePath(input$fapath),
                                         paired=input$read.paired,
                                         script=normalizePath(paste0(input$scrpath,'/events/run_featurecounts.sh')))
    config[['Task']][['generate_event_id']]=list(script=normalizePath(paste0(input$scrpath,'/events/generate_event_id.R')),
                                                 remove_chr=input$remove_chr,
                                                 gtf_path=normalizePath(input$GTFpath),
                                                 event_type=paste(input$event.types,collapse = ';')
                                                 )
    config[['Task']][['event_to_split_read_matrix']]=list(script=normalizePath(paste0(input$scrpath,'/events/event_to_split_read_matrix.r')),
                                                          ft_event_path=normalizePath(paste0(config$Basic$scrpath,'/data/gtex/')))
    config[['Task']][['reads_to_psi']]=list(script=normalizePath(paste0(input$scrpath,'/events/reads_to_psi.r')))
    config[['Task']][['velocity']]=list(output_path=normalizePath(paste0(config$Basic$work_path,'/velocity')))

    dir.create(config$Basic$work_path,recursive = T)    
    write(x = toJSON(x = config,pretty = T,auto_unbox = T),
          file = normalizePath(paste0(config$Basic$work_path,'/',config$DataSet,'.json')))
  })
}

