function combineDistance(inpath,outpath)
    [rbp,method,event_types]=parseInput(inpath);
    cors=struct();
    for i=1:length(event_types)
        type=event_types{i};
        disp([type,' ',datestr(datetime())])
        embedding=h5read(inpath,['/feature/',type]);
        psi=h5read(inpath,['/psi/',type]);
        if size(rbp,1)>=3
            %distance=comDist(rbp,psi,embedding,method);
            %%
            disp(['Calculate splicing regulation information ',datestr(datetime())]);
            correlation=single(corr(psi',rbp'));
            correlation(isnan(correlation))=0;
            clear 'psi';
            
            %disp(['Get pvalue ',datestr(datetime())]);
            [dist0,pvalue]=corr(correlation');
            clear 'dist0';
            pvalue(isnan(pvalue))=1;

            %disp(['Get alpha ',datestr(datetime())]);
            b=log(9999);
            a=-2*b;
            alpha=1./(1+exp(a*pvalue+b));
            clear 'pvalue';

            disp(['Get splicing regulation distance ',datestr(datetime())]);
            dist1=pdist(correlation,method);
            clear 'correlation';
            dist1=squareform(dist1);
            disp(['Get sequence distance ',datestr(datetime())]);
            dist2=pdist(embedding,method);
            clear 'embedding' ;
            dist2=squareform(dist2);
            
            disp(['Dist Normalize ',datestr(datetime())]);
            min1=min(dist1,[],'all');
            min2=min(dist2,[],'all');
            max1=max(dist1,[],'all');
            max2=max(dist1,[],'all');
            dist1=(dist1-min1)/(max1-min1);
            dist2=(dist2-min2)/(max2-min2);
            
            disp(['Combine Distance ',datestr(datetime())]);
            distance=(1-alpha).*dist1+alpha.*dist2;
            clear 'dist1' 'dist2';
        else
            disp(['Get sequence distance ',datestr(datetime())]);
            dist1=pdist(embedding,method);
            clear 'embedding' ;
            dist1=squareform(dist1);
            
            disp(['Dist Normalize ',datestr(datetime())]);
            min1=min(dist1,[],'all');
            max1=max(dist1,[],'all');
            dist1=(dist1-min1)/(max1-min1);
            
            disp(['Get Distance ',datestr(datetime())]);
            distance=dist1;
            clear 'dist1';
        end
        %%
        disp(['Save similarity ',type,' ',datestr(datetime())]);
        h5create(outpath,['/',type],size(distance));
        h5write(outpath,['/',type],distance);
    end
end

function [rbp,method,event_types]=parseInput(inpath)
    method=h5read(inpath,'/method');
    method=method{1};
    rbp=h5read(inpath,'/rbp');
    info=h5info(inpath);
    event_types={};
    for i=1:length(info.Groups(1).Datasets)
        event_types{i}=info.Groups(1).Datasets(i).Name;
    end
end
