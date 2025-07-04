showSelect=function(ff)
{
  if(class(input[[ff]])[1]!="list")
  {
    return(F)
  }

  if(.GlobalEnv$file.attr[ff,'type']=='shinyFilesButton')
  {
    path=parseFilePaths(volumns,input[[ff]])
    path=normalizePath(path$datapath)
  }
  else
  {
    path=parseDirPath(volumns,input[[ff]])
    path=normalizePath(path)
  }
  updateTextInput(session = session,inputId = gsub(pattern = "_",replacement = "",x = ff),value = path)
  # updateValue(gsub(pattern = "_",replacement = "",x = ff),env)
  return(T)
}

getDefaultPath=function(cmd)
{
  sys=Sys.info()['sysname']
  if(sys=='Windows')
  {
    msg=system(paste('where',cmd),intern = T,wait = T)
  }
  else
  {
    msg=system(paste('which',cmd),intern = T,wait = T)
  }

  if(!is.null(attr(msg,"status"))&&attr(msg,'status')==1)
  {
    return("")
  }
  return(normalizePath(msg))
}

getVolumes=function(exclude=NULL)
{
  osSystem <- Sys.info()["sysname"]
  if (osSystem == "Darwin") {
    volumes <- dir_ls("/Volumes")
    names(volumes) <- basename(volumes)
  }
  else if (osSystem == "Linux") {
    volumes <- c(Computer = "/")
    if (isTRUE(dir_exists("/media"))) {
      media <- dir_ls("/media")
      names(media) <- basename(media)
      volumes <- c(volumes, media)
    }
  }
  else if (osSystem == "Windows") {
    wmic <- paste0(Sys.getenv("SystemRoot"), "\\System32\\Wbem\\WMIC.exe")
    if (!file.exists(wmic)) {
      volumes_info <- system2("powershell", "$dvr=[System.IO.DriveInfo]::GetDrives();Write-Output $dvr.length $dvr.name $dvr.VolumeLabel;",
                              stdout = TRUE)
      num = as.integer(volumes_info[1])
      if (num == 0)
        return(NULL)
      mat <- matrix(volumes_info[-1], nrow = num, ncol = 2)
      mat[, 1] <- gsub(":\\\\$", ":/", mat[, 1])
      sel <- mat[, 2] == ""
      mat[sel, 2] <- mat[sel, 1]
      volumes <- mat[, 1]
      volNames <- mat[, 2]
      volNames <- paste0(volNames, " (", gsub(":/$", ":",
                                              volumes), ")")
    }
    else {
      volumes <- system(paste(wmic, "logicaldisk get Caption"),
                        intern = TRUE, ignore.stderr = TRUE)
      volumes <- sub(" *\\r$", "", volumes)
      keep <- !tolower(volumes) %in% c("caption", "")
      volumes <- volumes[keep]
      volNames <- system(paste(paste('cmd /c chcp 936 &&'),wmic, "/FAILFAST:1000 logicaldisk get VolumeName"),
                         intern = TRUE, ignore.stderr = TRUE)[-1]
      volNames = iconv(volNames,from = 'GBK',to = 'UTF8')
      volNames <- sub(" *\\r$", "", volNames)
      volNames <- volNames[keep]
      volNames <- paste0(volNames, ifelse(volNames == "",
                                          "", " "))
      volNames <- paste0(volNames, "(", volumes, ")")
    }
    names(volumes) <- volNames
    volumes <- gsub(":$", ":/", volumes)
  }
  else {
    stop("unsupported OS")
  }
  if (!is.null(exclude)) {
    volumes <- volumes[!names(volumes) %in% exclude]
  }
  volumes
}

initValue=function()
{
  browser()
  isolate({
    values=reactiveValuesToList(input)
  })
  return(values)
}

updateValue=function(id)
{
  isolate({
    new.value=input[[id]]
  })
  flag=identical(class(new.value),class(raw.values[[id]]))
  if(flag)
  {
    flag=identical(new.value,raw.values[[id]])
  }

  if(!flag)
  {
    raw.values[[id]]<<-new.value
  }
  printMSG(paste('updateValue:',id,flag))
  return(flag)
}

printMSG=function(cmd)
{
  time=Sys.time()
  cmd=paste(paste0('[',time,']'),cmd)
  print(cmd)
}

checkMsg=function(id,name,type,msg)
{
  flag=T
  msg=paste(msg,paste("checking",name,'...'))
  updateTextAreaInput(session = session,inputId = 'check_configure',label = '',value = msg)
  if(input[[id]]=='')
  {
    error.list=data.frame(id=id,msg=name,type=type,status='Error')
    flag=F
    msg=paste(msg,'Error!\r\n')
    updateTextAreaInput(session = session,inputId = 'check_configure',label = '',value = msg)
  }
  else
  {
    error.list=data.frame(id=id,msg=name,type=type,status='Success')
    msg=paste(msg,'Success!\r\n')
    updateTextAreaInput(session = session,inputId = 'check_configure',label = '',value = msg)
  }
  return(list(flag=flag,msg=msg,error=error.list))
}

validation=function(input,session)
{
  show_alert(title = 'Checking Configuration',
                 text = tags$span(id="check",
                                 textAreaInput(inputId = 'check_configure',label = '',value = "",width = '100%',height = '100%',resize = 'both',rows = 20),
                                 tags$h4(id='check_msg','')),
                 closeOnClickOutside = F,html = T,showCloseButton = F)
  flag=T
  msg=""

  # checklist=data.frame(id=c('dataset','cellinfofile','bampath','SCESCsrc','Rscriptpath','JAVApath','Samtoolspath','featurecountspath',
  #                           'MCRpath','rmatspath','MAJIQpath','IRFinderpath','workpath','ref_name','fapath','GTFpath','GFFpath',
  #                           'RBPpath','finetunerawpath','ref.pkg','phastpath'),
  #                      name=c('Dataset Name','Cell Info File','Bam Path','SCESC Source Directory','Rscript Path','JAVA Path','Samtools Path',
  #                             'FeatureCounts Path','MCR Path','rMats Path','MAJIQ Path','IRFinder Path','Work Path','Reference Name',
  #                             'Reference Genome','GTF File','GFF File','RBP File','Fine-tune Raw File','Reference Package','PhastCons Path'),
  #                      type=c('text','text','text','text','text','text','text','text','text','text','text','text',
  #                             'text','text','text','text','text','text','text','text','text'))

  checklist = data.frame(
    id = c(
      "dataset", "bampath", "condabinpath", "conda_envname", "pythonpath", "JAVApath", "Samtoolspath", "featurecountspath",
      "MCRpath", "rmatspath", "MAJIQ_env", "MAJIQlicensepath", "IRFinderpath",
      "STARpath", "workpath", "ref_name", "fapath", "GTFpath", "GFFpath",
      "RBPpath", "phastpath"
    ),
    name = c(
      "Dataset Name", "Bam Path", "Conda bin Path", "SCSES conda env", "Python Path", "JAVA Path", "Samtools Path",
      "FeatureCounts Path", "MCR Path", "rMats Path", "MAJIQ env",
      "MAJIQ license file", "IRFinder Path", "STAR Path", "Work Path", "Reference Name",
      "Reference Genome", "GTF File", "GFF File", "RBP File", "PhastCons Path"
    ),
    type = c(
      "text","text", "text", "text", "text", "text", "text", "text", "text", "text", "text",
      "text", "text", "text", "text", "text", "text", "text", "text", "text", "text"
    )
  )

  error_list=data.frame()
  for(i in seq(1,nrow(checklist)))
  {
    check.result=checkMsg(id = checklist$id[i],name = checklist$name[i],type=checklist$type[i],msg=msg)
    flag=flag&check.result$flag
    msg=check.result$msg
    error=check.result$error
    session$sendCustomMessage('checkStatus', as.list(error))
    error_list=rbind(error_list,error)
  }
  if("Error" %in% error_list$status){
    session$sendCustomMessage('checkResponse', list(status="fail"))
  }else{
    session$sendCustomMessage('checkResponse', list(status="success"))
  }
  # closeSweetAlert(session)
  return(flag)
}
