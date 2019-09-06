# sc_Matlab_toolbox
A toolbox for single cell sequencing data analysis (sc_Matlab_toolbox).
It covers the procedures of sc-RNA data analysis including 1), normalization for cells with a scaling factor (TMM baised w); 2,filteration for genes which have a lot of missed value in expression matrix; 3), identification of highly variable genes (HVGs); 4) clustering and visualization with t-SNE; and 5) identify the differential expression genes (DEGs) with three kinds of test (NB, Poisson, KW-test).
In the future, more will be added.
Requirement: Version of Matlab is more than 2018 to enshure having "tsne" function
Usage:
1. Retrieve data from https://www.ebi.ac.uk/gxa/sc/experiments?species=homo%20sapiens. You will get a text file in which there are three files. '.MTX' is matrix of which rows is gene index and cols is cell index. Another two files '.MTX_COLS' and '.MTX_ROWS' indicate the gene name corresponding the gene index and the cell ID to cell index, respectively.
2. Read the files you downloaded using "read_cell_atlas.m" by inputting the path wher you saved the three files. The function "read_cell_atlas.m" will generate vairable E of which "E.M " is matrix (row is gene and col is cell), E.cols=cell_name, and E.rows=gene_name. Then save the variable as one file (give name, suppose 'matrix_scRNAE_GEOD_83139.mat' at current path).
3. Call fcuntion "step_all.m" to run analsyis.
