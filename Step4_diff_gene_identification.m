clear;
clc; close all;
load ('S_Mscaled_matrix_scRNA.mat'); % S and M_scaled
load ('matrix_scRNA.mat'); % E.rows
% load ('S_Mscaled_matrix_scRNAE_GEOD_83139.mat');
% load('matrix_scRNAE_GEOD_83139.mat')
gene=E.rows;
Total_cluster=unique(S); %cluster number
cluster_num=length(Total_cluster);
meth_type=3;
isplot=1;
[Sort_M_scaled Sort_gene]=sort_gene_by_percent(M_scaled,gene);
M_scaled=Sort_M_scaled;gene=Sort_gene;


[DG]=my_Identify_diff_gene_one2others(M_scaled,S,gene,meth_type,isplot);
%[DG]=my_Identify_diff_gene_pairwise(M_scaled,S,gene,meth_type,isplot);


