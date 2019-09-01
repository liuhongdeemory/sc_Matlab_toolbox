function Y = my_tSNE(data_b,dis_type,Perplexity);
opt_type={'mahalanobis';'cosine';'chebychev';'euclidean'};
rng('default') % for reproducibility
if nargin<3
    Perplexity=25;
end
if nargin<2
    dis_type='cosine';
end
if isnan(find(strcmp(opt_type,dis_type)==1))
    dis_type='cosine';
    disp('Optinal distance type:')    
    disp(opt_type)
    return;
end

Y = tsne(data_b,'Distance',dis_type,'Perplexity',Perplexity);
