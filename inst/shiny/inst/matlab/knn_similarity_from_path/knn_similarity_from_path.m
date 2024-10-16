function knn_similarity_from_path(inpath,outpath)
    [similar,type,paras]=parseInput(inpath);
    ka=ceil(paras.k/3);
    sigma=mink(similar,ka+1,2);
    sigma=sigma(:,ka+1);
    similar=-1*(similar./sigma).^2;
    similar=exp(similar);
    
    similar=(similar+similar')/2;
    knn_thresh=maxk(similar,paras.k+1,2);
    knn_thresh=knn_thresh(:,paras.k+1);
    
    similar(similar<knn_thresh)=0;
    similar=sparse(similar);
    
    similar=similar./sum(similar,2);
    similar(find(isnan(similar)))=0;
    m1=similar;
    delta=Inf;
    step=0;
    while(delta>paras.decay)
        m2=paras.alpha*m1*similar+(1-paras.alpha)*similar;
        sse=sum(sum((m2-m1).^2));
        ssr=full(m2)-mean(m1,2);
        ssr=ssr.^2;
        ssr=sum(sum(ssr));
        delta=sse/(sse+ssr);
        disp(delta);
        m1=m2;
        step=step+1;
        if(step>50)
            break
        end
    end
    clear m2
    knn_thresh=maxk(m1,paras.k+1,2);
    knn_thresh=knn_thresh(:,paras.k+1);
    m1(m1<knn_thresh)=0;
    m1=full(m1);
    similar=m1./sum(m1,2);
    similar(find(isnan(similar)))=0;
    similar=sparse(similar);
    save(outpath,'similar');
end
function [data,type,paras]=parseInput(inpath)
    type=h5read(inpath,'/type');
    type=type{1};
    dist_path=h5read(inpath,'/dist_path');
    dist_path=dist_path{1};
    data=h5read(dist_path,['/',type]);
    paras=struct();
    paras.k=h5read(inpath,'/paras/k');
    paras.alpha=h5read(inpath,'/paras/alpha');
    paras.decay=h5read(inpath,'/paras/decay');
end