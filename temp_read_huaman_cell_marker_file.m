clear all; clc;
% CM= readtable('human_cell_marker_M.xls');
% marker.tissueType=CM{:,2};
% marker.UberonOntologyID=CM{:,3};
% marker.cancerType=CM{:,4};
% marker.cellType=CM{:,5};
% marker.cellName=CM{:,6};
% marker.CellOntologyID=CM{:,7};
% marker.cellMarker=CM{:,8};
% marker.geneSymbol=CM{:,9};
% marker.geneID=CM{:,10};
% marker.proteinName=CM{:,11};
% marker.proteinID=CM{:,12};
% marker.markerResource=CM{:,13};
% marker.PMID=CM{:,14};
load ('hg_cell_marker.mat');
for i=1:length(marker.geneSymbol)
    s=marker.geneSymbol{i};
    S=regexp(s, ',\s', 'split');
    marker.eachSymbol(i).marker=S;
end
