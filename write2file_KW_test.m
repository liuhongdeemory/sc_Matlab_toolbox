function write2file_KW_test(gene,chi_s,P,FDR,file_in)
[gene,chi_s,P,FDR]=sort_by_pValue_KW_test(gene,chi_s,P,FDR);
gene=EnsID2geneSymbol(gene);
if length(P)>=1
    if ~isnan(file_in)
        % sorting by FDR or P
        fid = fopen(file_in, 'w');

        disp(strcat('Total diff gene:',num2str(length(gene))))
        fprintf(fid, '%s\t %s\t %s\t %s\t\n',...
            'Gene Symbol','Chi-square', '-log10(P)', 'FDR');
        % write title
        for i=1:length(gene)
            g=gene{i};
            chi=chi_s(i);
            p=P(i);
            fdr=FDR(i);
            fprintf(fid, '%s\t %d\t %d\t %d\t\n',...
                g, chi, p,  fdr);
        end
        
        fclose(fid);
        disp('done')
    else
        disp('File path is empty.')
    end
    disp('File path is empty.')
else
    disp('DEG are Not found.')
    return;
end
