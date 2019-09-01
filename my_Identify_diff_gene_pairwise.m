function [DG]=my_Identify_diff_gene_pairwise(M_scaled,S,gene,meth_type,isplot)
if nargin<3
    disp('Variables M_scaled,S and gene are needed')
    return
elseif nargin==3
    meth_type=1;
    isplot=0;
elseif nargin==4
    isplot=0;
end
gene=EnsID2geneSymbol(gene);
Total_cluster=unique(S); %cluster number
cluster_num=length(Total_cluster);
DG=[]; % record the Diff gene
disp('Identify differential expression gene:')
switch meth_type
    case 1
        disp('method=Negative Bionomial Model')
        k=1;
        for i=1:cluster_num-1
            for j=i+1:cluster_num
                mean_i=mean(M_scaled(:,find(S==i)),2);
                mean_j=mean(M_scaled(:,find(S==j)),2);
                FC=log2(mean_i./mean_j);% FC
                tLocal = nbintest(M_scaled(:,find(S==i)),M_scaled(:,find(S==j)),'VarianceLink','LocalRegression');
                p=tLocal.pValue;
                nlg_p=-log10(p);
                [FDR] = mafdr(p,'Showplot', false);% FDR estimation
                % write into files
                file_in=strcat('NB Model for-',num2str(i),'-VS-',num2str(j),'.txt');
                
                write2file_NB_P(gene,mean_i,mean_j,FC,nlg_p,FDR,file_in);
                % disp(gene_index)
                if isplot==2
                    plotVarianceLink(tLocal)
                    title(file_in,'fontsize',14)
                end
                
                % write into a struct array
                [s_gene,s_mean_i,s_mean_j,s_FC,s_nlg_p,s_FDR]=sort_by_pValue_NB(gene,mean_i,mean_j,FC,nlg_p,FDR);
                DG(k).header=file_in;
                DG(k).gene=s_gene';
                temp=[s_mean_i s_mean_j s_FC s_nlg_p s_FDR];
                DG(k).stat=temp;
                k=k+1;
            end
        end
        
    case 2
        disp('method=Possion Model')
        header={};gene_index=[];num_group=[];k=1;
        for i=1:cluster_num-1
            for j=i+1:cluster_num
                mean_i=mean(M_scaled(:,find(S==i)),2);
                mean_j=mean(M_scaled(:,find(S==j)),2);
                FC=log2(mean_i./mean_j);% FC
                tLocal = nbintest(M_scaled(:,find(S==i)),M_scaled(:,find(S==j)),'VarianceLink','Identity');
                p=tLocal.pValue;
                nlg_p=-log10(p);
                [FDR] = mafdr(p,'Showplot', false);% FDR estimation
                % write into files
                file_in=strcat('Possion Model for-',num2str(i),'-VS-',num2str(j),'.txt');
                write2file_NB_P(gene,mean_i,mean_j,FC,nlg_p,FDR,file_in);
                if isplot==2
                    plotVarianceLink(tLocal)
                    title(file_in,'fontsize',14)
                end
                
                % write into a struct array
                [s_gene,s_mean_i,s_mean_j,s_FC,s_nlg_p,s_FDR]=sort_by_pValue_NB(gene,mean_i,mean_j,FC,nlg_p,FDR);
                DG(k).header=file_in;
                DG(k).gene=s_gene';
                temp=[s_mean_i s_mean_j s_FC s_nlg_p s_FDR];
                DG(k).stat=temp;
                k=k+1;
            end
        end
        
    case 3
        disp('method=Kruskal-Wallis test.')
        header={};gene_index=[];num_group=[];k=1;
        M_scaled=10000*log10(M_scaled+1);
        for i=1:length(gene)
            [p(i) tbl] = kruskalwallis(M_scaled(i,:),S,'off');
            nlg_p(i)=-log10(p(i));
            chi_s(i)= tbl{2,5};
        end
        [FDR] = mafdr(p,'Showplot', false);% FDR estimation
        file_in=strcat('KW-test.txt');
        write2file_KW_test(gene,chi_s,nlg_p,FDR,file_in);
        
        % write into a struct array
        [s_gene,s_chi_s,s_nlg_p,s_FDR]=sort_by_pValue_KW_test(gene,chi_s,nlg_p,FDR);
        DG(k).header=file_in;
        DG(k).gene=s_gene';
        temp=[s_chi_s' s_nlg_p' s_FDR'];
        DG(k).stat=temp;
        k=k+1;
        
    otherwise
        disp('Err method.')
        return;
end
disp('Files are saved in current path.')

%%%%% sort unique DG gene
if isplot==1
    S_sort=[];
    for i=1:cluster_num
        temp= find(S==i);
        S_sort=[S_sort;temp];
    end
    G_sort=[];
    for i=1:length(DG)
        temp=find_gene_idx(gene,DG(i).gene);
        temp=temp(~isnan(temp));
        G_sort=[G_sort;temp];
    end
    my_heatmap(log2(M_scaled(G_sort,S_sort)+1));
    %%% setting figure
    colormap('hot')
    set(gca,'ylim',[-5 length(G_sort)]);
    set(gca,'xlim',[1 length(S_sort)+15]);
    hold on;
    c=hsv(cluster_num);
    x(1)=1;
    for i=1:cluster_num
        x(i+1)=length(find(S<=i));
        plot([x(i) x(i+1)],[-2 -2],'linewidth',5,'color',c(i,:));
        posi_x=round((x(i)+x(i+1))/2);
        text(posi_x,-3,strcat('Cluster-',num2str(i)),'fontsize',14,'rotation',45)
    end 
    G=[];
    y(1)=1;
    for i=1:length(DG)-1
        temp=find_gene_idx(gene,DG(i).gene); % only marker the first gene
        temp=temp(~isnan(temp));
        G=[G;temp];        
        y(i+1)=length(G)+1;        
    end
    for i=1:length(y)
        text(length(S_sort)+1,y(i),DG(i).gene{1},'fontsize',10)
    end
    axis off
end
