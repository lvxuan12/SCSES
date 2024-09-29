function imputation1(datapath,outpath)
    disp(datetime)
    [similar,data_type,similar_type,parameter,data]=parseInput(datapath);
    disp(datetime)
    fields=fieldnames(data);
    data=struct2cell(data);
    if strcmp(data_type,'RC')
        for i=1:size(data,1)
            data{i}=run(similar,similar_type,data{i}',parameter.decay)';
        end
    else
        for i=1:size(data,1)
            data{i}=run(similar,similar_type,data{i}',parameter.decay)';
        end
    end
    disp(datetime)
    data=cell2struct(data,fields,1);
    save(outpath,'data','-v7.3');
    disp(datetime)
end

function [similar,data_type,similar_type,parameter,data]=parseInput(datapath)
    similar=h5read(datapath,'/similar');
    similar=sparse(similar);
    data_type=h5read(datapath,'/data_type');
    data_type=data_type{1};
    similar_type=h5read(datapath,'/similar_type');
    similar_type=similar_type{1};
    parameter=struct();
    parameter.decay=h5read(datapath,'/parameter/decay');
    dataInfo=h5info(datapath,'/data');
    data=struct();
    fields=strings;
    for i=1:size(dataInfo.Datasets,1)
        field=dataInfo.Datasets(i).Name;
        fields{i}=field;
        value=h5read(datapath,strcat('/data/',field));
        data=setfield(data,field,value);
    end
end

function new=run(similar,similar_type,data,decay)
    if(strcmp(similar_type,'event'))
        data=data';
    end    
    error=Inf;
    new=data;
    step=1;
    while(error>decay)
        current=new;
        new=similar*current;
%        error=new-current;
%        error=mean(error.*error,'all');
        sse=sum(sum((new-current).^2));
        ssr=sum(sum((new-mean(current,2)).^2));
        error=sse/(sse+ssr);
        disp(error);
        
        step=step+1;
        if(step>50)
            break
        end
    end
    if(strcmp(similar_type,'event'))
        new=new';
    end
end
