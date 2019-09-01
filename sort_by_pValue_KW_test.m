function [s_gene,s_chi_s,s_P,s_FDR]=sort_by_pValue_KW_test(gene,chi_s,P,FDR);
top_num=10;
[i_]=find(P>=6 & FDR<=0.01);
if ~isempty(i_)
    %disp('There are genes which satisfy the limit.')
    s_gene=gene(i_);
    s_chi_s=chi_s(i_);
    s_P=P(i_);
    s_FDR=FDR(i_);
    
    [temp, j_]=sort(s_FDR,'ascend');
    s_gene=s_gene(j_);
    s_chi_s=s_chi_s(j_);
    s_P=s_P(j_);
    s_FDR=s_FDR(j_);
else
    disp('There is no gene taouchs the limit; shown is top 10 by FDR.')
    [temp, j_]=sort(FDR,'ascend');
    s_gene=gene(j_(1:top_num));
    s_chi_s=chi_s(j_(1:top_num));
    s_P=P(j_(1:top_num));
    s_FDR=FDR(j_(1:top_num));
end

