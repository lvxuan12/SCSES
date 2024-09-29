Shiny.addCustomMessageHandler("checkStatus",function(msg) {
    if(msg.status=='Success')
    {
      if(msg.type=='text')
      {
        $("#"+msg.id).css('border-color','#5cb85c')
      }
    }
    else if(msg.status=='Error')
    {
      if(msg.type=='text')
      {
        $("#"+msg.id).css('border-color','#d9534f')
      }
    }
});
