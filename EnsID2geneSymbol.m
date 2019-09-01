function geneSymbol=EnsID2geneSymbol(EnsID)
load('hg_ENS_symbol.mat');
geneSymbol={};

for i=1:length(EnsID)
    i_=find(strcmp(hg.ENS,EnsID(i)));
    if ~isnan(i_)
        geneSymbol{i}=hg.symbol{i_(1)};
    else
        i_=find(strcmp(hg.symbol,EnsID(i)));
        if ~isnan(i_)
            geneSymbol{i}=hg.symbol{i_(1)};
        else
            geneSymbol{i}=EnsID{i};
        end        
    end
end
geneSymbol=geneSymbol';

clear hg;