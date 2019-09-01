function [Sort_M_scaled Sort_gene]=sort_gene_by_percent(M_scaled,gene,cut)
per=[];
if nargin <3
cut=0.8;
end

[gene_num cell_num]=size(M_scaled);
for i=1:gene_num
    per(i)=sum(M_scaled(i,:)==0)/cell_num;
end

[i_]=find(per>cut);
Sort_M_scaled=M_scaled(i_,:);
Sort_gene=gene(i_);
disp(strcat('Number of genes expressed in more than', num2str(cut),'% cells:'))
disp(length(Sort_gene))
