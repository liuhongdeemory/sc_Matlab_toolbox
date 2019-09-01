function write2file_NB_P(gene,mean_i,mean_j,FC,P,FDR,file_in)
[s_gene,s_mean_i,s_mean_j,s_FC,s_P,s_FDR]=sort_by_pValue_NB(gene,mean_i,mean_j,FC,P,FDR);
s_gene=EnsID2geneSymbol(s_gene);
if length(P)>=1
    if ~isnan(file_in)
        fid = fopen(file_in, 'w');
        
        disp(strcat('Total diff gene:',num2str(length(s_gene))))
        fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t\n',...
            'Gene Symbol', 'Mean-i','Mean-j','log2(FC)', '-log10(P)', 'FDR');
        % write title
        for i=1:length(s_gene)
            g=s_gene{i};
            mi=s_mean_i(i);
            mj=s_mean_j(i);
            fc=s_FC(i);
            p=s_P(i);
            fdr=s_FDR(i);
            fprintf(fid, '%s\t %d\t %d\t %d\t %d\t %d\t\n',...
                g, mi, mj, fc, p,  fdr);
        end
        
        fclose(fid);
        disp('done')
    else
        disp('File path is empty.')
    end
else
    disp('DEG are Not found.')
    return;
end
