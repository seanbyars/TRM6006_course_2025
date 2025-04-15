# scRNAseq analysis outputs for Assignment 1b
- This provides an example of the main outputs required from your single cell analysis.
- Once you generate these outputs using your chosen single cell dataset (e.g. four available, see 'single_cell_datasets_for_assignment_use.r'), these can go in your report.
- these outputs were generated using the demo dataset 'tutorial_example_seurat_combined.Rds' and the R code in 'Assignment 1b instructions.r'

# List of outputs
1. single cell preprocessing
- violin plot of cell QC metrics before filtering [see example_output1.png]
- violin plot of cell QC metrics after filtering [see example_output2.png]
- UMAP showing cell clustering before integration [see example_output3.png]
- UMAP showing cell clustering after integration [see example_output4.png]
- UMAP showing your predicted cell types [see example_output5.png]
2. DGE cluster biomarker analysis
- violin plot of 2-3 significant genes of interest from one cell type [see example_output6.png]
- feature plot of these same 2-3 genes [see example_output7.png]
- top 20 significant genes in a table (for one cell type) [see example_output8.png]
3. DGE control-treat analysis
- violin plots of 2-3 genes of interest, showing differences in cases vs control for one cell type [see example_output9.png]
- feature plots of these same 2-3 genes, showing differences in cases vs controls for one cell type [see example_output10.png]
- top 20 genes in a table (for one cell type) [see example_output11.png]
4. gene ontology enrichment analysis
- at a minimum, a table of the top 10-20 significant/interesting GO terms (include GO term, number of genes in GO term, number of intersecting genes, fold enrichment, FDr/adjusted p-value) [see example_output12.png]

