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
  volumns=c(Home=path_home(),getVolumes())
  env=environment()

  raw.values=isolate(reactiveValuesToList(input))

  for(ff in file.attr$id)
  {
    if(file.attr[ff,'type']=='shinyFilesButton')
    {
      shinyFileChoose(input = input, id = ff, roots=volumns,session = session,defaultRoot = 'Home')
    }
    else
    {
      shinyDirChoose(input = input, id = ff, roots=volumns,session = session,defaultRoot = 'Home')
    }
  }

  updateTextInput(session = session, inputId = gsub(pattern = "_", replacement = "", x = "python_path"), value = getDefaultPath(cmd = "python"))
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'JAVA_path'),value = getDefaultPath(cmd = 'java'))
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'Samtools_path'),value = getDefaultPath(cmd = 'samtools'))
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'featurecounts_path'),value = getDefaultPath(cmd = 'featureCounts'))
  updateTextInput(session = session, inputId = gsub(pattern = "_", replacement = "", x = "IRFinder_path"), value = getDefaultPath(cmd = "IRFinder"))
  updateTextInput(session = session, inputId = gsub(pattern = "_", replacement = "", x = "STAR_path"), value = getDefaultPath(cmd = "STAR"))
  updateTextInput(session = session, inputId = gsub(pattern = "_", replacement = "", x = "rmats_path"), value = getDefaultPath(cmd = "rmats.py"))
# updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = 'Matlab_path'),value = getDefaultPath(cmd = 'Matlab_path'))

  # textinput action----
  # observeEvent(c(input$cellinfo_file,input$cellinfofile),{
  #     if(showSelect('cellinfo_file'))
  #     {
  #       isolate({
  #         cellinfofile=input$cellinfofile
  #       })
  #       updateValue("cellinfo_file")
  #       if(!updateValue("cellinfofile"))
  #       {
  #         if(cellinfofile!="")
  #         {
  #           if(grepl(pattern = "\\.csv$",x = cellinfofile))
  #           {
  #             sep=','
  #           }
  #           else
  #           {
  #             sep="\t"
  #           }
  #           tryCatch({
  #             printMSG(paste('Read New File',cellinfofile))
  #             data=read.table(cellinfofile,header = T,sep = sep,check.names = F)
  #             columns=colnames(data)
  #             names(columns)=columns
  #             updateRadioGroupButtons(session = session,inputId = 'cellinfo.column',choices = columns,selected = columns[1])
  #           },error=function(e){
  #             print(e)
  #             updateRadioGroupButtons(session = session,inputId = 'cellinfo.column',choices = 'NA',selected = 'NA')
  #           })
  #
  #         }
  #       }
  #     }
  # })
  observeEvent(c(input$bam_path,input$bampath),{
    showSelect('bam_path')
    updateValue("bam_path")
    updateValue("bampath")

  })
  observeEvent(c(input$python_path,input$pythonpath),{
    showSelect('python_path')
    updateValue("python_path")
    updateValue("pythonpath")
  })

  observeEvent(c(input$JAVA_path,input$JAVApath),{
    showSelect('JAVA_path')
    updateValue("JAVA_path")
    updateValue("JAVApath")
  })
  observeEvent(c(input$MCR_path,input$MCRpath),{
    showSelect('MCR_path')
    updateValue("MCR_path")
    updateValue("MCRpath")
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
  # observeEvent(c(input$Matlab_path,input$Matlabpath),{
  #   showSelect('Matlab_path')
  #   updateValue("Matlab_path")
  #   updateValue("Matlabpath")
  # })
  observeEvent(c(input$rmats_path,input$rmatspath),{
    showSelect('rmats_path')
    updateValue("rmats_path")
    updateValue("rmatspath")
  })
  # observeEvent(c(input$MAJIQ_path,input$MAJIQpath),{
  #   showSelect('MAJIQ_path')
  #   updateValue("MAJIQ_path")
  #   updateValue("MAJIQpath")
  # })
  observeEvent(c(input$MAJIQ_license_path,input$MAJIQlicensepath),{
    showSelect('MAJIQ_license_path')
    updateValue("MAJIQ_license_path")
    updateValue("MAJIQlicensepath")
  })
  # observeEvent(c(input$VOILA_path,input$VOILApath),{
  #   showSelect('VOILA_path')
  #   updateValue("VOILA_path")
  #   updateValue("VOILApath")
  # })
  observeEvent(c(input$IRFinder_path,input$IRFinderpath),{
    showSelect('IRFinder_path')
    updateValue("IRFinder_path")
    updateValue("IRFinderpath")
  })
  observeEvent(c(input$STAR_path,input$STARpath),{
    showSelect('STAR_path')
    updateValue("STAR_path")
    updateValue("STARpath")
  })
  # observeEvent(c(input$BRIE_path,input$BRIEpath),{
  #   showSelect('BRIE_path')
  #   updateValue("BRIE_path")
  #   updateValue("BRIEpath")
  # })
  # observeEvent(c(input$outrigger_path,input$outriggerpath),{
  #   showSelect('outrigger_path')
  #   updateValue("outrigger_path")
  #   updateValue("outriggerpath")
  # })
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
  # observeEvent(c(input$finetune_raw_path,input$finetunerawpath),{
  #   showSelect('finetune_raw_path')
  #   updateValue("finetune_raw_path")
  #   updateValue("finetunerawpath")
  # })
  observeEvent(c(input$phast_path,input$phastpath),{
    showSelect('phast_path')
    updateValue("phast_path")
    updateValue("phastpath")
  })
  # observeEvent(c(input$brie_ref,input$brieref),{
  #   showSelect('brie_ref')
  #   updateValue("brie_ref")
  #   updateValue("brieref")
  # })
  # observeEvent(c(input$brie_feature,input$briefeature),{
  #   showSelect('brie_feature')
  #   updateValue("brie_feature")
  #   updateValue("briefeature")
  # })
  # observeEvent(c(input$outrigger_size,input$outriggersize),{
  #   showSelect('outrigger_size')
  #   updateValue("outrigger_size")
  #   updateValue("outriggersize")
  # })

  # create configuration action----
  observeEvent(input$create_configure,{
    config=list()
    flag=validation(input,session)
    if(!flag)
    {
      return()
    }
    config[['DataSet']]=input$dataset
    config[['Basic']][['sequence']]=input$read.umi
    config[['Basic']][['readlength']]=input$readlength
    config[['Basic']][['paired']]=input$read.paired
    config[['Basic']][['bam_path']]=normalizePath(input$bampath)
    config[['Basic']][['work_path']]=normalizePath(input$workpath)
    config[['Basic']][['core']]=input$max.core
    config[["Basic"]][["conda_envname"]] = input$conda_envname
    config[["Basic"]][["java_path"]] = normalizePath(input$JAVApath)
    config[["Basic"]][["python_path"]] = normalizePath(input$pythonpath)
    config[['Basic']][['mcr_path']]=normalizePath(input$MCRpath)
    config[["Basic"]][["STAR_path"]] = normalizePath(input$STARpath)
    config[["Basic"]][["samtools_path"]] = normalizePath(input$Samtoolspath)
    config[["Basic"]][["featureCounts_path"]] = normalizePath(input$featurecountspath)
    config[["Basic"]][["rMATS_path"]] = normalizePath(input$rmatspath)
    config[["Basic"]][["IRFinder_path"]] = normalizePath(input$IRFinderpath)
    config[['Basic']][['MAJIQ_env']]=input$MAJIQ_env
    config[["Basic"]][["refgenome"]][["gtf_path"]] = normalizePath(input$GTFpath)
    config[["Basic"]][["refgenome"]][["gff_path"]] = normalizePath(input$GFFpath)
    config[["Basic"]][["refgenome"]][["ref_path"]] = normalizePath(input$fapath)
    config[["Basic"]][["refgenome"]][["genome_name"]] = input$ref_name

    config[['Basic']][['filter_merged_bam']]=list(ExonToIntronReads=input$exon.intron.read,junctionReads=input$junction.read)
    config[["Basic"]][["filter_events_sc"]] = list(
      minCell = input$minCell, minRC = input$minRC,
      min.percentCells.gene = input$exp.dropout.ratio,
      min.percentCells.event = input$psi.dropout.ratio,
      min.nFeatures = input$cell.gene.exp, min.nCount = input$cell.read,
      max.percentMT = input$cell.mt.pct,
      filter.mt = input$filter.mt, filter.rp = input$filter.rp
    )

    config[["Task"]][["event"]][["remove_chr"]] = input$remove_chr
    config[["Task"]][["event"]][["event_type"]] = input$event.types
    config[["Task"]][["event"]][["majiq_license_file"]] = input$MAJIQ_license_path

    config[["Task"]][["impute"]][["rbp"]] = normalizePath(input$RBPpath)
    config[["Task"]][["impute"]][["cell_similarity_data"]] = paste(input$cell.similar, collapse = ";")
    config[["Task"]][["impute"]][["feature_num"]] = input$feature_num

    config[['Task']][['impute']][['event_features']]=list(phast_path=normalizePath(input$phastpath),
                                                          chr_prefix=input$chr.prefix)

    for(type in input$event.types)
    {
      config[['Task']][['impute']][['event_features']][['AE']][[type]]=list(epoch=input[[paste0(type,'.epoch')]],
                                                 embedding=input[[paste0(type,'.embedding')]],
                                                 layer=paste0('[',input[[paste0(type,'.layer')]],']')
                                                )
    }
    config[["Task"]][["impute"]][["KNN"]][["cell"]] = list(
      distance_method = input$similar.method,
      kmax = input$cell.kman,
      kmin = input$cell.kmin,
      alpha = input$cell.alpha,
      decay = input$cell.decay
    )
    config[["Task"]][["impute"]][["KNN"]][["event"]] = list(
      k = input$event.k,
      alpha = input$event.alpha,
      decay = input$event.decay
    )
    config[["Task"]][["impute"]][["decay_impute"]] = input$all.decay

    dir.create(config$Basic$work_path,recursive = T)
    write(x = toJSON(x = config,pretty = T,auto_unbox = T),
          file = normalizePath(paste0(config$Basic$work_path,'/',config$DataSet,'.json')))
  })
}

