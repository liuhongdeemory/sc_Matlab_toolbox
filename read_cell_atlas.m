function [E]=read_cell_atlas(path_)
F_=dir(path_);

for J=3:length(F_)
    if strcmpi(F_(J).name(end-3:end),'.MTX')==1 & ~isnan(findstr(F_(J).name,'normalised'))
        f_mtx=strcat(path_,'\',F_(J).name);
    end
    if strcmpi(F_(J).name(end-8:end),'.MTX_COLS')==1 & ~isnan(findstr(F_(J).name,'normalised'))
        f_mtx_cols=strcat(path_,'\',F_(J).name);
    end
    if strcmpi(F_(J).name(end-8:end),'.MTX_ROWS')==1 & ~isnan(findstr(F_(J).name,'normalised'))
        f_mtx_rows=strcat(path_,'\',F_(J).name);
    end
end



disp(f_mtx)
disp(f_mtx_cols)
disp(f_mtx_cols)

%%%%%%% read file
disp('open maxtrix data ....')
fid1 = fopen(f_mtx);
i=0;
while feof(fid1)~=1 %%& i<=1999
    tline = fgetl(fid1);
    
    S=regexp(tline, '\s', 'split');
    if i>=3
        gene_id(i-2)=str2num(S{1});% get gene id
        cell_id(i-2)=str2num(S{2});% get cell id
        count_(i-2)=str2num(S{3}); % expression
    end
    i=i+1;
end
fclose(fid1);
S = sparse(gene_id, cell_id,count_);
disp('done.')
% construct matrix
M_gene=unique(gene_id);M_cell=unique(cell_id);
%行和列的标记
row_n=length(M_gene);  col_n=length(M_cell);
M=zeros(row_n,col_n);
disp('construct matrix ....')
M=full(S);
disp('done')
%%%%%%% read rows
fid1 = fopen(f_mtx_rows);
i=1;
while feof(fid1)~=1 %& i<=1999
    tline = fgetl(fid1);
    S=regexp(tline, '\s', 'split');
    gene_name(i)=S(1);% get gene id
    i=i+1;
end
fclose(fid1);
gene_name=gene_name(M_gene);

%%%%%%% read cols
fid1 = fopen(f_mtx_cols);
i=1;
while feof(fid1)~=1 %& i<=1999
    tline = fgetl(fid1);    
    S=regexp(tline, '\s', 'split');
    cell_name(i)=S(1);% get gene id
    i=i+1;
end
fclose(fid1);
cell_name=cell_name(M_cell);

E.M=M;
E.cols=cell_name;
E.rows=gene_name;
save('matrix_scRNAE_GEOD_83139.mat','E')
disp('matrix_scRNA.mat is saved.')
%E:\PBMC_single_cell\E-GEOD-83139-normalised-files

