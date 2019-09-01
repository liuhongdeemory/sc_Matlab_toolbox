function [M_scaled,w]=scale_normalization(M,cutoff,isplot)
if nargin<2
    cutoff=0.1;
end
if nargin<3
    isplot=0;
end

[gene_num cell_num]=size(M);
% 10 times for test
m_val=[];
for J=1:10
    r_idx=randperm(cell_num);
    cut_site=round(cell_num/2);
    ref=M(:,[r_idx(1:cut_site)]);
    test=M(:,[cut_site+1:cell_num]);
    for i=1:gene_num
        av_ref=mean(ref(i,:));
        av_test=mean(test(i,:));
        m_val(i,J)=log2(av_test+1)-log2(av_ref+1);
    end
end
mm_val=mean(m_val');
i_=find(abs(mm_val)>=0.01 & abs(mm_val)<=cutoff);
disp(strcat('Number of less-changed genes:',num2str(length(i_))));
if ~isnan(i_)
    a=mean(M(i_,:));
    w=a(1)*(1./a); % this scale factor
    
    for i=1:gene_num
    M_scaled(i,:)=M(i,:).*w;
    end
else
    disp('please set a proper cutoff.')
    return
end
%%ploting
if isplot==1 
    subplot(1,2,1)
    histogram(mm_val,100,'normalization','probability','facecolor',[1,1,1],'edgecolor',[0 0 0])
    hold on
    plot([-0.01 -0.01],[0 0.4],'b:','linewidth',0.5);
    plot([0.01 0.01],[0 0.4],'b:','linewidth',0.5);
    plot([-cutoff -cutoff],[0 0.4],'r:','linewidth',0.5);
    plot([cutoff cutoff],[0 0.4],'r:','linewidth',0.5);
    xlabel('log2(test/ref)','fontsize',14)
    ylabel('Gene count','fontsize',14)
    set(gca,'fontsize',14)
    set(gca,'ylim',[-0.01 0.4],'xlim',[min(mm_val) max(mm_val)])
    
    subplot(2,2,2)
    counts=M(i_,1:10);
    maboxplot(log10(counts),'title','Raw Read Count','orientation','horizontal')
    ylabel('sample','fontsize',14)
    xlabel('counts(log10)','fontsize',14)    
    subplot(2,2,4)
    normCounts=M_scaled(i_,1:10);
    maboxplot(log10(normCounts),'title','Normalized Read Count','orientation','horizontal')
    ylabel('sample','fontsize',14)
    xlabel('counts(log10)','fontsize',14)   
end

