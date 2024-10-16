function scses(datapath,outpath)
    disp(datetime)
    [similar,parameter,data]=parseInput(datapath,'PSI','cell');
    data=struct2cell(data);
    psi_psi_cell=run(similar,'cell',data{1}',parameter.decay)';
    
    [similar,parameter,data]=parseInput(datapath,'RC','cell');
    fields=fieldnames(data);
    data=struct2cell(data);
    for i=1:size(data,1)
        data{i}=run(similar,'cell',data{i}',parameter.decay)';
    end
    idx_in=find(contains(fields,'retention')==1);
    idx_ex=find(contains(fields,'exclusion')==1);
    in_reads=zeros(size(data{1},1),size(data{1},2));
    ex_reads=zeros(size(data{1},1),size(data{1},2));
    for i=1:size(idx_in,2)
        disp(fields{idx_in(i)});
        in_reads=in_reads+data{idx_in(i)};
    end
    for i=1:size(idx_ex,2)
        disp(fields{idx_ex(i)});
        ex_reads=ex_reads+data{idx_ex(i)};
    end
    mean_reads=(in_reads+ex_reads)/(size(idx_in,2)+size(idx_ex,2));
    in_reads=in_reads/size(idx_in,2);
    ex_reads=ex_reads/size(idx_ex,2);
    psi_rc_cell=in_reads./(in_reads+ex_reads);
    psi_rc_cell(isnan(psi_rc_cell))=0;

    [row,col,v_in] = find(in_reads);
    min_reads_in=prctile(v_in,10);
    [row,col,v_ex] = find(ex_reads);
    min_reads_ex=prctile(v_ex,10);

    flag=in_reads>min_reads_in;
    psi_rc_cell=psi_rc_cell.*flag;
    psi_rc_cell(find((in_reads>min_reads_in)&(ex_reads<min_reads_ex)))=1;
    
    [similar,parameter,data]=parseInput(datapath,'PSI','event');
    disp(size(psi_rc_cell,1));
    if(size(psi_rc_cell,1)==1)
        data=struct('v1_psi',psi_psi_cell,'v1_rc',psi_rc_cell,'v2_psi',psi_rc_cell);
    else
        psi_psi_event=run(similar,'event',psi_rc_cell',parameter.decay)';
        data=struct('v1_psi',psi_psi_cell,'v1_rc',psi_rc_cell,'v2_psi',psi_psi_event);
    end
    save(outpath,'data','-v7.3');
    disp(datetime)
end

function [similar,parameter,data]=parseInput(datapath,data_type,similar_type)
    similar=h5read(datapath,['/similar/',similar_type]);
    similar=sparse(similar);
    parameter=struct();
    parameter.decay=h5read(datapath,'/parameter/decay');
    dataInfo=h5info(datapath,['/data/',data_type]);
    data=struct();
    fields=strings;
    for i=1:size(dataInfo.Datasets,1)
        field=dataInfo.Datasets(i).Name;
        fields{i}=field;
        value=h5read(datapath,strcat('/data/',data_type,'/',field));
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
        if(step>50 || size(current,2)==1)
            break
        end
    end
    if(strcmp(similar_type,'event'))
        new=new';
    end
end
