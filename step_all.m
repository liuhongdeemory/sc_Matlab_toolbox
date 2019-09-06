clear all; clc; close all;
disp('load data.')
load ('matrix_scRNAE_GEOD_83139.mat');
% load ('matrix_scRNA.mat');
X=E.M;gene=E.rows; % original data
isplot=1;

disp('Scaling for each cell.')
cutoff=0.1;
figure1=figure('color','w');
[M,w]=scale_normalization(X,cutoff,isplot); % M is the scaled X; and w is factor;
M=X;

disp('filtering genes.')
cut=0.4;
[Sort_M Sort_gene]=sort_gene_by_percent(M,gene,cut); % you can skip this function
% if the step is skipped, the inputted data for "my_Find_HVGs.m" should be "M" and "gene".

disp('Find highly variable genes.')
figure1=figure('color','w');
[X_HVGs,HVGs]=my_Find_HVGs(Sort_M,Sort_gene,isplot);

% X_HVGs=Sort_M;HVGs=Sort_gene;

disp('tSNE.')
isplot=1;
dis_type={'chebychev';'cosine';'euclidean';'mahalanobis'};
Perplexity=40; % para
Y = my_tSNE(X_HVGs',dis_type{2},Perplexity);
% if you would not like to use martix of highly variable genes, you can use M to replace X_HVGs.


disp('Identify subgourps for tSNE (clustering).')
figure1=figure('color','w');
cluster_num=[];isplot=1;
opt_type=[1:4];%opt_type_full={'genlouvain';'louvain';'kmeans';'linkage'};
[S Q]=my_Cluster_algorithm(Y(:,1),Y(:,2),opt_type(1),cluster_num,isplot);
% S is community (cluster) ID for each cell.
disp(strcat('Q=',num2str(Q)));
disp(strcat('Num of clusters=',num2str(length(unique(S)))));

disp('identifying differentially expressed genes and write them in to files');
meth_type=1;isplot=1;
figure1=figure('color','w');
[DG]=my_Identify_diff_gene_one2others(X_HVGs,S,HVGs,meth_type,isplot);
