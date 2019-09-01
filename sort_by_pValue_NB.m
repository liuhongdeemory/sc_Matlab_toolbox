function [s_gene,s_mean_i,s_mean_j,s_FC,s_P,s_FDR]=sort_by_pValue_NB(gene,mean_i,mean_j,FC,P,FDR);
% sorting by P and report top 10 genes
top_num=10;
[i_]=find(P>=6 & FC>=2 & FDR<=0.01);
if ~isempty(i_)
    %disp('There are genes which satisfy the limit.')
    
    s_gene=gene(i_);
    s_mean_i=mean_i(i_);
    s_mean_j=mean_j(i_);
    s_FC=FC(i_);
    s_P=P(i_);
    s_FDR=FDR(i_);
    
    [temp, j_]=sort(s_FDR,'ascend');
    s_gene=s_gene(j_);
    s_mean_i=s_mean_i(j_);
    s_mean_j=s_mean_j(j_);
    s_FC=s_FC(j_);
    s_P=s_P(j_);
    s_FDR=s_FDR(j_);
else
    disp('There is no gene touchs the limit; shown is top 10 by FDR.')
    
    [temp, j_]=sort(FDR,'ascend');
    s_gene=gene(j_(1:top_num));
    s_mean_i=mean_i(j_(1:top_num));
    s_mean_j=mean_j(j_(1:top_num));
    s_FC=FC(j_(1:top_num));
    s_P=P(j_(1:top_num));
    s_FDR=FDR(j_(1:top_num));
end



