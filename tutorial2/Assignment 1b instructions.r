# ASSIGNMENT 1b DATA and ANALYSIS INSTRUCTIONS


# DEMO DATA & CODE
# the example below provides an R code walkthrough single cell analysis with the Lupus data (2 x Lupus vs 2 x control samples) from the Friday tutorial
# This starts with loading an single cell Rds file into R, then performs some standard single cell processing (e.g. QC filtering, normalisation) and analysis (cluster biomarker analysis, DGE analysis)
# If you would like to download this data and run through each step below on your laptop, the data is available at https://tinyurl.com/5y9buyce and will download a file called 'tutorial_example_seurat_combined.Rds'


# DATA TO USE FOR YOUR ASSIGNMENT
# You will need to choose a different dataset for your single cell analysis
# I have provided FOUR different single cell datasets - you can choose one of these for your assignment, or you are also welcome to choose your own from GEO, SRA or another source, but you will still need to include the main outputs below in your report
# Datasets for use - see the document 'single_cell_datasets_for_assignment_use.r' for download instructions and further information on each dataset


# OUTPUTS TO INCLUDE FOR YOUR ASSIGNMENT
# There are 4 main analysis parts (1. single cell preprocessing, 2. cluster biomarker analysis, 3. DGE control-treat analysis, 4. GO analysis) to complete
# I've listed below the main outputs needed to go into your assignment, they are:

# 1. single cell preprocessing
#	- violin plot of cell QC metrics before filtering [OUPUT FOR REPORT 1]
#	- violin plot of cell QC metrics after filtering [OUPUT FOR REPORT 2]
#	- UMAP showing cell clustering before integration [OUPUT FOR REPORT 3]
#	- UMAP showing cell clustering after integration [OUPUT FOR REPORT 4]
#	- UMAP showing your predicted cell types [OUPUT FOR REPORT 5]
# 2. DGE cluster biomarker analysis
#	- violin plot of 2-3 significant genes of interest from one cell type [OUPUT FOR REPORT 6]
#	- feature plot of these same 2-3 genes [OUPUT FOR REPORT 7]
#	- top 20 significant genes in a table (for one cell type) [OUPUT FOR REPORT 8]
# 3. DGE control-treat analysis
#	- violin plots of 2-3 genes of interest, showing differences in cases vs control for one cell type [OUPUT FOR REPORT 9]
#	- feature plots of these same 2-3 genes, showing differences in cases vs controls for one cell type [OUPUT FOR REPORT 10]
#	- top 20 genes in a table (for one cell type) [OUPUT FOR REPORT 11]
# 4. gene ontology enrichment analysis
#	- at a minimum, a table of the top 10-20 significant/interesting GO terms (include GO term, number of genes in GO term, number of intersecting genes, fold enrichment, FDr/adjusted p-value) [OUPUT FOR REPORT 12]
#	- you could also use some visualisations, below I've provided some links to different online tools specialised in GO visualisation

# In addition to this, please provide your R code - you can copy and paste this into a word document, a text or R file, or include within an R markdown report if you are using Rstudio




# R CODE USED WITH DEMO DATASET + description of outputs needed



# Load your single cell dataset



library(Seurat)

setwd("~/Desktop/demo_data")
seurat_combined=LoadSeuratRds("tutorial_example_seurat_combined.Rds")

# have a look at the info in the Seurat object - here we can see 41,902 cells, with 32,833 genes detected, in four samples (which are presently in different layers)
seurat_combined
# An object of class Seurat 
# 32833 features across 41902 samples within 4 assays 
# Active assay: RNA (32738 features, 0 variable features)
#  4 layers present: counts.1, counts.2, counts.3, counts.4
#  3 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3





# 1. single cell preprocessing



# what does the quality of the cells look like BEFORE filtering?
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01) # [OUPUT FOR REPORT 1]

# here I'm applying some QC filters - note these need to be adjusted specifically for your data
# the ones that I have used below may need to be changed for your data
# have a look at the bioinformatics platform workshop here - https://monashbioinformaticsplatform.github.io/scRNAseq_Workshop_ABACBS_2024/qc.html
# they used a slightly different filtering, e.g. seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt < 5)
dim(seurat_combined)
seurat_combined <- subset(seurat_combined, 
      subset = nFeature_RNA <= 2500 & 
      nFeature_RNA > 200 &
      percent.mt <= 15)
dim(seurat_combined)

# what does the quality of the cells look like AFTER filtering?
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01) # [OUPUT FOR REPORT 2]

# Normalize data - adjusts for sequencing depth per cell
seurat_combined <- NormalizeData(seurat_combined)

# Identify variable features - finds biological signal in the data - used in downstream analysis, e.g. PCA
seurat_combined <- FindVariableFeatures(seurat_combined)

# Scale the data - gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
seurat_combined <- ScaleData(seurat_combined)

# Run PCA - dimensionality reduction for visualisation of cell-cell similarity
seurat_combined <- RunPCA(seurat_combined, npcs = 30)

# FindNeighbors() function is used to compute a nearest-neighbor graph of cells based on a dimensionality reduction (typically PCA). This graph is then used for clustering (via FindClusters()) and visualization.
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:30, reduction = "pca")

# FindClusters() used to identify clusters of cells based on the neighbor graph created with FindNeighbors(). It performs community detection using algorithms like Louvain, SLM, or Leiden.
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5, cluster.name = "unintegrated_clusters")

# RunUMAP() computes a UMAP (Uniform Manifold Approximation and Projection) dimensional reduction for visualization of cells in 2D or 3D space, usually after PCA and clustering.
seurat_combined <- RunUMAP(seurat_combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# have a look at the UMAP - are there clusters specific to each sample?
# The cells across the four samples are not yet integrated
# Integration is needed as it combines datasets from different experiments, batches, or conditions to correct for technical variation and align shared biological signals. This enables meaningful comparisons and joint analysis across cells from multiple sources.
DimPlot(seurat_combined, reduction = "umap.unintegrated") # this is a UMAP coloured by cluster - the number of clusters here is determined by the cluster resolution selected above - different clusters typically represent different cell types.
DimPlot(seurat_combined, reduction = "umap.unintegrated", group.by="Sample") # [OUPUT FOR REPORT 3] colour the cells by what sample they come from
DimPlot(seurat_combined, reduction = "umap.unintegrated", split.by="Sample", ncol=2) # we can also split the UMAP into different samples

# Now we need to run the integration - harmony is one memory-efficient way, but there are several others, e.g. see https://satijalab.org/seurat/articles/seurat5_integration
# If you have not yet installed the R harmony package, see https://github.com/immunogenomics/harmony for install instructions
library(harmony)

seurat_combined <- RunHarmony(
  object = seurat_combined,
  group.by.vars = "Sample",  # Specify the 'Sample' column to integrate by - you may need to change this if your Sample column has a different name.
  reduction.use = "pca",        # Use PCA as the reduction method
  dims = 1:30               # Number of dimensions to use
)

# Find neighbors and clusters
seurat_combined <- FindNeighbors(seurat_combined, reduction = "harmony", dims = 1:30)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5) # we're just using a resolution of 0.5. 

# Update embeddings to use Harmony
seurat_combined <- RunUMAP(seurat_combined, reduction = "harmony", dims = 1:30)

# Now we can have a look at cell clustering again with the UMAP, with the integrated data this time
DimPlot(seurat_combined, reduction = "umap")
DimPlot(seurat_combined, reduction = "umap", group.by = "Sample") # [OUPUT FOR REPORT 4]
DimPlot(seurat_combined, reduction = "umap", split.by = "Sample", ncol=2)
# We can also look at the clustering by predicted cell type - you may need to change this depending on what predicted cell types you have in your dataset
# In this example below, because I previously predicted celltypes using Azimuth's PBMC reference, I have three different cell type predictions, level 1 (broad) through to level 3 (specific) cell types
DimPlot(seurat_combined, reduction = "umap", group.by = "predicted.celltype.l1") # [OUPUT FOR REPORT 5]


# minimum outputs to include in your report from '1. single cell preprocessing'
- violin plot of cell QC metrics before filtering [OUPUT FOR REPORT 1]
- violin plot of cell QC metrics after filtering [OUPUT FOR REPORT 2]
- UMAP showing cell clustering before integration [OUPUT FOR REPORT 3]
- UMAP showing cell clustering after integration [OUPUT FOR REPORT 4]
- UMAP showing your predicted cell types [OUPUT FOR REPORT 5]





# 2. DGE cluster biomarker analysis



# From here, we can look at what genes are differentially expressed in each cell type, relative to all other cell types
# This analysis is also called a cluster biomarker analysis, as it tells us how gene expression differs in each cell type
# Have a look at the FindAllMarkers() Seurat function first
?FindAllMarkers
# set the cell type variable that we want to run the cluster biomarker analysis on
Idents(seurat_combined) <- "predicted.celltype.l1"
# If you have layers in your Seurat object, you may also need to run this function
seurat_combined=JoinLayers(seurat_combined)
# now run the DGE analysis
de_result_wilcox <- FindAllMarkers(seurat_combined)
# have a look at the results
head(de_result_wilcox)
dim(de_result_wilcox)

# for this analysis, choose some differentially expressed genes that you're interested in and plot them by cell type (e.g. Violin plot) and also look at their expression on the UMAP
# In this example, I'm interested in a couple of genes (KLRF1,IGFBP7) as I could see they were highly significantly differentially expressed in NK cells
de_result_wilcox=de_result_wilcox[rev(order(de_result_wilcox$avg_log2FC)),] # ordering the dataset by fold change..
de_result_wilcox=de_result_wilcox[(order(de_result_wilcox$p_val_adj)),] # then FDR value
de_result_wilcox_NK=de_result_wilcox[c(which(de_result_wilcox$cluster == "NK")),] # and now just pulling out results for NK cells
de_result_wilcox_NK[1:10,]
#          p_val avg_log2FC pct.1 pct.2 p_val_adj cluster    gene
#SH2D1B.3      0   6.810671 0.255 0.004         0      NK  SH2D1B
#BOK           0   6.611722 0.066 0.001         0      NK     BOK
#LGALS9B       0   6.605866 0.074 0.001         0      NK LGALS9B
#KLRF1.4       0   6.428588 0.796 0.022         0      NK   KLRF1 - KLRF1, much more highly expressed in NK cells
#RAMP1.1       0   6.279516 0.111 0.002         0      NK   RAMP1
#AKR1C3.4      0   5.451664 0.302 0.012         0      NK  AKR1C3
#NCAM1.1       0   5.283878 0.092 0.003         0      NK   NCAM1
#IGFBP7.4      0   5.261569 0.425 0.029         0      NK  IGFBP7 - IGFBP7, also more highly expressed in NK cells
#IL18RAP.4     0   5.170031 0.160 0.006         0      NK IL18RAP
#B3GNT7.1      0   5.073871 0.106 0.005         0      NK  B3GNT7

VlnPlot(seurat_combined, features = c("KLRF1","IGFBP7")) # [OUPUT FOR REPORT 6] violin plot showing expression of these genes - it will plot this for predicted.celltype.l1 cell types as we set this above
FeaturePlot(seurat_combined, reduction = "umap", features = c("KLRF1","IGFBP7")) # [OUPUT FOR REPORT 7]
DimPlot(seurat_combined, reduction = "umap",label=TRUE)  # Visualize by cell type (predicted.celltype.l1)

# [OUPUT FOR REPORT 8] write out the top 20 genes from this analysis, for your chosen cell type
write.csv(de_result_wilcox_NK[1:20,],"~/Desktop/1_cluster_biomarker_analysis_NK_cells_top20genes.csv", row.names=FALSE)

# outputs to include in your report from this analysis - I've labelled these above as [OUPUT FOR REPORT]
- violin plot of 2-3 genes of interest from above [OUPUT FOR REPORT 6]
- feature plot of those same genes [OUPUT FOR REPORT 7]
- top 20 genes in a table [OUPUT FOR REPORT 8]





# 3. DGE control-treat analysis



# The next thing we can look at is if there are genes that are significantly different in certain cell type between cases (Lupus) vs. controls
# we first need to define control and lupus samples in our data - we do this by adding an additional column to our Seurat metadata, which can then be used in the differential gene expression analysis
META=seurat_combined@meta.data # pulling out the Seurat metadata into a separate dataframe
META$CASE_CONTROL="CASE" # create new column
META$CASE_CONTROL[c(which(META$Sample == "cHD1"))]="CONTROL" # now code those that are controls in our dataset, there are just two control samples (no Lupus)
META$CASE_CONTROL[c(which(META$Sample == "cHD2"))]="CONTROL"
table(META$Sample, META$CASE_CONTROL) # check that this worked properly, there should be cells for both cases and controls
#          CASE CONTROL
#  cHD1       0   13276
#  cHD2       0    8706
#  cSLE10 13290       0
#  cSLE27  5665       0
META=META[,c("Sample","CASE_CONTROL")] # now just selecting the columns I want to add back into the Seurat metadata
seurat_combined = AddMetaData(seurat_combined, META) # we use the AddMetaData() function to add this
head(seurat_combined) # you should now see this new column in your Seurat metadata

# anno level 1 DGE - case vs control

# Because we are running a differential gene expression analysis between cases and controls, for each cell type, we need to creat one more column in the Seurat metadata that defines cases and controls for each celltype. I'm just going to add this column directly.
seurat_combined$LEVEL1.treat <- paste(seurat_combined$predicted.celltype.l1, seurat_combined$CASE_CONTROL, sep = "_")
head(seurat_combined) # you should now see an extra column in the Seurat metadata
table(seurat_combined$LEVEL1.treat) # lets have a look at the new column - so there are 2398 Lupus B cells, 3290 control B cells as an example
#         B_CASE       B_CONTROL      CD4 T_CASE   CD4 T_CONTROL      CD8 T_CASE   CD8 T_CONTROL         DC_CASE      DC_CONTROL       Mono_CASE    Mono_CONTROL         NK_CASE      NK_CONTROL    other T_CASE other T_CONTROL 
#           2275            3239            6925            7892            2899            5962              66             229            5457            2660             787            1239             386             705 
#     other_CASE   other_CONTROL 
#            160              56 
# Because we are using this new column to run the DGE analysis, we point to this with Idents()
Idents(seurat_combined) <- "LEVEL1.treat"
LEVELS=names(table(seurat_combined$predicted.celltype.l1)) # this is just a list of the cell types. I'm just going to run this analysis with the level 1 predicted cell types.
for(i in 1:length(LEVELS)) { # now using a loop in R, to cycle through each cell type, testing for gene expression differences between Cases vs Controls
  print(i)
  SLE=paste(LEVELS[i],"_CASE",sep=""); SLE
  CONTROL=paste(LEVELS[i],"_CONTROL",sep=""); CONTROL
  print(SLE); print(CONTROL)
  CLUSTERX_control_treat <- FindMarkers(seurat_combined, ident.1 = SLE, ident.2 = CONTROL, verbose = FALSE)
  CLUSTERX_control_treat$cluster=LEVELS[i]
  CLUSTERX_control_treat$gene=rownames(CLUSTERX_control_treat)
  if(i == 1) {DGE_analysis_case_vs_control=CLUSTERX_control_treat} else {DGE_analysis_case_vs_control=rbind(DGE_analysis_case_vs_control,CLUSTERX_control_treat)}
}
# have a look at the results output
head(DGE_analysis_case_vs_control)
dim(DGE_analysis_case_vs_control)

# this has run a DGE analysis between cases (Lupus) vs control, for each cell type
# here I'm just going to focus on NK cells again

DGE_analysis_case_vs_control_NK=DGE_analysis_case_vs_control[c(which(DGE_analysis_case_vs_control$cluster == "NK")),]
dim(DGE_analysis_case_vs_control_NK)
DGE_analysis_case_vs_control_NK[1:10,] # just looking at the top 10
#                p_val avg_log2FC pct.1 pct.2     p_val_adj cluster   gene
#RPS265  1.099592e-145 -1.3775726 0.841 0.978 3.599843e-141      NK  RPS26
#IFI44L5 6.106004e-125  5.4150543 0.450 0.028 1.998984e-120      NK IFI44L
#LGALS15 3.640068e-104  2.0935464 0.809 0.395  1.191685e-99      NK LGALS1
#GNLY5   1.247841e-102 -0.7253820 0.991 1.000  4.085180e-98      NK   GNLY
#HBB5     1.982303e-98 -0.2194414 0.371 0.021  6.489665e-94      NK    HBB
#IFI65    2.175982e-88  3.2718599 0.549 0.169  7.123731e-84      NK   IFI6
#ISG205   3.313398e-88  1.7985003 0.732 0.326  1.084740e-83      NK  ISG20
#IFITM35  3.164522e-84  1.9628176 0.663 0.249  1.036001e-79      NK IFITM3
#IL325    5.814965e-83  1.6499854 0.775 0.361  1.903703e-78      NK   IL32
#MX15     4.786599e-81  3.8987866 0.347 0.034  1.567037e-76      NK    MX1
# I can see that genes like IFI44L and LGALS1 are significantly upregulated in Lupus NK cells, relative to control NK cells
# I'm going to do some visualisation with these genes
VlnPlot(seurat_combined, features = c("IFI44L"), split.by = "CASE_CONTROL", group.by = "predicted.celltype.l1") # [OUPUT FOR REPORT 9] violin plot, showing the expression of these genes in cases vs controls
VlnPlot(seurat_combined, features = c("LGALS1"), split.by = "CASE_CONTROL", group.by = "predicted.celltype.l1") # [OUPUT FOR REPORT 9]
FeaturePlot(seurat_combined, reduction = "umap", features = c("IFI44L"), split.by="CASE_CONTROL") # [OUPUT FOR REPORT 10] UMAP plot showing how the expression of these genes differs between cases vs controls
FeaturePlot(seurat_combined, reduction = "umap", features = c("LGALS1"), split.by="CASE_CONTROL") # [OUPUT FOR REPORT 10]

# just checking the expression of these two genes in all other cell types that we tested
# IFI44L also shows differences between cases vs controls for other cell types
DGE_analysis_case_vs_control_NK[c(which(DGE_analysis_case_vs_control_NK$gene == "IFI44L"))] # OUT[c(which(OUT$gene == "IFI44L")),]
#                p_val avg_log2FC pct.1 pct.2     p_val_adj cluster   gene
#IFI44L   0.000000e+00   4.173040 0.492 0.060  0.000000e+00       B IFI44L
#IFI44L1  0.000000e+00   3.604638 0.401 0.074  0.000000e+00   CD4 T IFI44L
#IFI44L2  0.000000e+00   4.667818 0.404 0.046  0.000000e+00   CD8 T IFI44L
#IFI44L3  5.476775e-33   3.669531 0.743 0.203  1.792987e-28      DC IFI44L
#IFI44L4 4.515968e-258   3.083097 0.492 0.124 1.478438e-253    Mono IFI44L
#IFI44L5 5.124208e-135   5.413970 0.472 0.029 1.677563e-130      NK IFI44L
#IFI44L6  6.004627e-02   1.402663 0.208 0.108  1.000000e+00   other IFI44L
#IFI44L7  3.738435e-26   3.320008 0.289 0.063  1.223889e-21 other T IFI44L

DGE_analysis_case_vs_control_NK[c(which(DGE_analysis_case_vs_control_NK$gene == "LGALS1"))] # OUT[c(which(OUT$gene == "LGALS1")),]
#                p_val avg_log2FC pct.1 pct.2     p_val_adj cluster   gene
#LGALS1  1.510542e-108  1.9096663 0.398 0.150 4.945214e-104       B LGALS1
#LGALS11  1.462926e-97  1.3487108 0.306 0.172  4.789328e-93   CD4 T LGALS1
#LGALS12  9.427228e-05  0.1750262 0.440 0.380  1.000000e+00   CD8 T LGALS1
#LGALS13  5.952374e-06  0.7060091 0.936 0.929  1.948688e-01      DC LGALS1
#LGALS14  0.000000e+00  0.8628346 0.993 0.977  0.000000e+00    Mono LGALS1
#LGALS15 1.988357e-109  2.0771377 0.818 0.398 6.509483e-105      NK LGALS1
#LGALS16  2.792555e-02 -0.1782118 0.332 0.492  1.000000e+00   other LGALS1
#LGALS17  2.111961e-27  2.1754550 0.430 0.154  6.914137e-23 other T LGALS1

# [OUPUT FOR REPORT 11] write out the top 20 genes from this analysis, for your chosen cell type
write.csv(DGE_analysis_case_vs_control_NK[1:20,],"~/Desktop/2_case_control_analysis_NK_cells_top20genes.csv", row.names=FALSE)

# write out a list of genes for the gene ontology analysis below
write.csv(DGE_analysis_case_vs_control_NK[1:500,], "~/Desktop/2_case_control_analysis_NK_cells_top500genes_forGOanalysis.csv", row.names=FALSE) # top 500 for geneontology.org tool
write.csv(DGE_analysis_case_vs_control_NK, "~/Desktop/2_case_control_analysis_NK_cells_allgenes_forGOrilla_GOanalysis.csv", row.names=FALSE) # writing out all genes if you want to use the GOrilla gene ontology tool

# outputs to include in your report from this analysis - I've labelled these above as [OUPUT FOR REPORT]
- violin plots of 2-3 genes of interest, showing differences in cases vs control [OUPUT FOR REPORT 9]
- feature plots of those same genes, showing differences in cases vs controls [OUPUT FOR REPORT 10]
- top 20 genes in a table [OUPUT FOR REPORT 11]





# 4. gene ontology enrichment analysis



# We are now going to look at what biological pathways are related to our set of differentially expressed genes
# For this analysis, choose one of your cell types from 2. above and write this out as a csv file - this can be used to run a gene ontology analysis of your choice
# I'm choosing NK cells again, my question coul be something like are biological pathways in Lupus NK cells perturbed, relative to control NK cells

# You can run gene ontology enrichment analysis a number of ways (in R using gene ontology enrichment tools, or using a web tool such as geneontology.org, GOrilla or others)
# the examples I gave in the Friday tutorial were to either:
# 1) choose the top few hundred (3-500) most differentially expressed genes and paste them into the gene ontology tool (https://geneontology.org/)
# 2) take the entire list of genes (ranked by p value) and paste them into the GOrilla gene ontology tool (https://cbl-gorilla.cs.technion.ac.il/)


# outputs to include in your report from this analysis
- for this analysis, you could present the significant/interesting GO pathways a number of different ways, for example:
- at a minimum, a table of the top 10-20 significant/interesting pathways (include GO term, number of genes in GO term, number of intersecting genes, fold enrichment, FDR/adjusted p-value) [OUPUT FOR REPORT 12]
- you could also use some visualisations, below I've provided some links to different online tools specialised in GO visualisation
https://bioinformatics.sdstate.edu/go/
https://geneontology.org/docs/tools-overview/
https://wego.genomics.cn/
http://revigo.irb.hr/
