file.attr=data.frame()
createTextFile=function(id,label,type='shinyFilesButton',...)
{
  ui=list(
    div(class="col-sm-2",style="padding-right:0",
        textInput(inputId = gsub(pattern = '_',replacement = "",x = id),label = label,...),
    ),
    div(class='col-sm-1',style='padding-left:0',
        div(class='form-group shiny-input-container',
            tags$label(class="control-label","placehold","for"=id,style="visibility: hidden"),
            get(type)(id = id,label = 'Browser',title = 'System Files',
                      multiple = F,buttonType = 'default',style = 'display:block;')
        )
    )
  )
  
  .GlobalEnv$file.attr=rbind(.GlobalEnv$file.attr,data.frame(id=id,type=type))
  rownames(.GlobalEnv$file.attr)=.GlobalEnv$file.attr$id
  return(ui)
}
