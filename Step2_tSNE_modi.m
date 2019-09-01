clear;
clc; close all;
% %%% perform tSNE
load ('matrix_scRNAE_GEOD_83139.mat');
%load ('matrix_scRNA.mat');
 cutoff=0.1;isplot=1;
[M_scaled,w]=scale_normalization(E.M,cutoff,isplot);


dis_type={'chebychev';'cosine';'euclidean';'mahalanobis'};Perplexity=25; % para
Y = my_tSNE(E.M',dis_type{2},Perplexity);
YY = my_tSNE(M_scaled',dis_type{2},Perplexity);
% %%% cluster the Y
opt_type_full={'genlouvain';'louvain';'kmeans';'linkage'};

x=Y(:,1);y= Y(:,2);
xx=YY(:,1);yy= YY(:,2);
opt_type=[1:4];
cluster_num=6;
isplot=1;
figure('color',[1 1 1])
for i=1:4
subplot(2,2,i)
[S Q(i)]=my_Cluster_algorithm(x,y,opt_type(i),cluster_num,isplot);

end
figure('color',[1 1 1])
for i=1:1
subplot(2,2,i)
[S Q(i)]=my_Cluster_algorithm(xx,yy,opt_type(i),cluster_num,isplot);
end


