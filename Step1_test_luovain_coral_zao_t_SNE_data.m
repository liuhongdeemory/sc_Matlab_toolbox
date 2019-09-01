clear all; clc; close all;
load ('coral_zao_tSNE_data.mat'); % variables x and y
X=[x y];

% G = graph(A~=0);
subplot(2,2,1)
idx = kmeans(X,3);
gscatter(x,y,idx)

subplot(2,2,2)
Z = linkage(X,'ward');
c = cluster(Z,'Maxclust',3);
gscatter(x,y,c)

subplot(2,2,3)
[A]=constuct_dis_Matrix(x,y);

A=A<mean(mean(A))/3;
k = full(sum(A));
twom = sum(k);

% B = @(v) ((A(:,v) - k'*k(v)/twom)); % original
 B = @(v) ((A(:,v) - k'*k(v)/twom));
%B = A - k'*k/twom;
[S,Q] = iterated_genlouvain(B);

gscatter(x,y,S)
% [my_C]=do_sub_graph(G,S);
% plot(G,'MarkerSize',7, 'NodeColor',my_C,'EdgeColor','w','NodeLabel',S,'Layout','force');

subplot(2,2,4)
[A]=constuct_dis_Matrix(x,y);
% [COMTY ending] = cluster_jl(A);
[COMTY ending]=cluster_jl_orientT(A);
L=COMTY.COM{:,1};
gscatter(x,y,L)


