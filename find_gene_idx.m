function idx=find_gene_idx(all_EnsID,EnsID)

idx=[];
for i=1:length(EnsID)
    i_=find(strcmp(all_EnsID,EnsID(i)));
    if ~isnan(i_)
    idx(i)=i_(1);
    else
    idx(i)=nan;    
    end
end
idx=idx';
