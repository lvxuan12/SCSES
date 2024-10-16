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
        $("#"+msg.id).css('border-color','#FF410D').css('border-width','3px')
        $("#check_configure").css('border-color','#FF410D').css('border-width','3px')
      }
    }
});

Shiny.addCustomMessageHandler("checkResponse",function(msg) {
    if(msg.status=='success')
    {
        $("#check_msg").html("Create successfully").css('color','#28a745')
    }
    else
    {
      $("#check_msg").html("Essential content missing").css('color','#dc3545')
    }
});
