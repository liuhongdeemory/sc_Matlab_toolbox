function [S Q]=my_Cluster_algorithm(x,y,meth_type,cluster_num,isplot)
opt_type=[1:4];
opt_type_full={'genlouvain';'louvain';'kmeans';'linkage'};
if nargin==1
    [a,b]=szie(x);
    if b==2
        x1=x(:,1);
        x2=x(:,2);
    else
        disp('Error: data should have two col')
        return;
    end
    meth_type=1;
elseif nargin==2
    meth_type=1;
    x1=x;x2=y;
elseif nargin>=3
    if isnan(find(opt_type,meth_type)==1)
        disp('Error: Optinal distance type:')
        disp(opt_type)
        return;
    end
    x1=x;x2=y;
end
if isnan(isplot)
    isplot==0;
end


d=3;
switch meth_type
    case 1
        [A]=constuct_dis_Matrix(x1,x2);
        A=A<mean(mean(A))/d;
        k = full(sum(A));
        twom = sum(k);
        B = @(v) ((A(:,v) - k'*k(v)/twom));
        %B = A - k'*k/twom;
        [S,Q] = iterated_genlouvain(B);
        Q = Q/twom;
    case 2
        [A]=constuct_dis_Matrix(x,y);
        A=A<mean(mean(A))/d;
        % [COMTY ending] = cluster_jl(A);
        [COMTY ending]=cluster_jl_orientT(A);
        [Q i_]=max(COMTY.MOD);
        S=COMTY.COM{:,i_};
    case 3
        S=kmeans([x1 x2],cluster_num);
        Q=nan;
    case 4
        Z = linkage([x1 x2],'ward');
        S = cluster(Z,'Maxclust',cluster_num);
        Q=nan;
    otherwise
        disp('Error: Optinal distance type:')
        disp(opt_type)
        return;
end

if isplot==1
    gscatter(x,y,S)
    xlabel('tSNE-1','fontsize',14)
    ylabel('tSNE-2','fontsize',14)
    title(opt_type_full{meth_type},'fontsize',14)
    set(gca,'fontsize',14)
end


