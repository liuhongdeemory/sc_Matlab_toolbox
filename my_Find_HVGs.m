function [X_HVGs,HVGs]=my_Find_HVGs(X,gene_,isplot);
% X=E.M;gene_=E.rows;%cell_=E.cols;
if nargin<3
    isplot=0;
end

u=nanmean(X,2);
cv=nanstd(X,0,2)./u;
m=size(X,1);
ydata=log10(cv);
xdata=u(~isnan(ydata));
lgcv=ydata(~isnan(ydata));
lgm=log10(xdata);
[lgm i_]=sort(lgm,'ascend');
lgcv=lgcv(i_);gene_=gene_(i_);
xSeq=linspace(min(lgm),max(lgm),100);
ySeq=zeros(length(xSeq),1);
% cal median
ySeq(1)=median(lgcv<=xSeq(1));
flag=zeros(1,length(lgcv));
idx=[];
for i=2:length(xSeq)
    i_=(lgm>xSeq(i-1) & lgm<=xSeq(i));
    cv_temp=lgcv(i_);
    ySeq(i)=median(cv_temp);
    j_=find(lgm>xSeq(i-1) & lgm<=xSeq(i));
    if ~isempty(j_)
        for v=1:length(j_)
            if lgcv(j_(v))>1.1*ySeq(i);%%%%%%%%%% Pls modify the 2
                idx=[idx j_(v)];
            end
        end
        
    end
end
HVGs=gene_(idx);
X_HVGs=X(idx,:);
disp(strcat('Total genes =',num2str(length(gene_))));
disp(strcat('Number of HVGs =',num2str(length(HVGs))));
%%% plot
if isplot==1
%     fig=figure('color','w');
    scatter(lgm,lgcv,'k.');
    hold on
    plot(xSeq,ySeq,'b-','linewidth',2)
    i_=find(flag==1);
    scatter(lgm(idx),lgcv(idx),'r.')
    xlabel('log10(mean)','fontsize',14);
    ylabel('log10(CV)','fontsize',14);
    title('Highly variable genes (HVGs)','fontsize',14)
    legend({'All genes';'Median line of CV';'HVGs'});
    set(gca,'fontsize',14,'color',[1 1 1])
end