# ------ study data 1
- This is single cell data from the same study I used for the Friday tutorial, however this is a different selection of control/Lupus samples
Study title: A single cell approach to map cellular subsets involved in Systemic Lupus Erythematosus (SLE) heterogeneity
Paper: https://www.nature.com/articles/s41590-020-0743-0
original GEO source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135779
	GSM4029915 	cHD3 [JB17018] - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4029915
	GSM4029922 	cHD4 [JB18069] - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4029922
	GSM4029910 	cSLE14 [JB17021] - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4029910
	GSM4029933 	cSLE31 [JB18080] - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4029933

# you can download the Rds data here with these four samples - this data has predicted cell types included
https://tinyurl.com/4fabwrbj

# read in the data in with LoadSeuratRds
seurat_combined = LoadSeuratRds("study1_seurat_combined.Rds") # 24037 cells across 4 samples, two main groups (Lupus vs control)



# ------ study data 2
Study title: Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19
Paper: https://pubmed.ncbi.nlm.nih.gov/32651212/
original source 1: https://cellxgene.cziscience.com/collections/4f889ffc-d4bc-4748-905b-8eb9db47a2ed
original source 2: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689

# you can download the Rds data here, I've selected a reduced dataset of COVID-19 vs normal samples - this data has predicted cell types included
https://tinyurl.com/5a6fja4n

# read in the data in with LoadSeuratRds
seurat_combined = LoadSeuratRds("study2_seurat_combined.Rds") # 49053 cells across 15 samples, two main groups (COVID-19 vs normal)



# ------ study data 3
Study title: Selective advantage of mutant stem cells in human clonal hematopoiesis is associated with attenuated response to inflammation and aging
Paper: https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(24)00207-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1934590924002078%3Fshowall%3Dtrue
Original sources:
	https://portal.sds.ox.ac.uk/articles/dataset/TARGET-seq_transcriptome_raw_counts/23576379?file=43521621
	https://portal.sds.ox.ac.uk/articles/dataset/TARGET-seq_metadata/23576262?file=43521618
	https://portal.sds.ox.ac.uk/articles/dataset/TARGET-seq_genotyping_data/23576421?file=43521585

# you can download the Rds data here - this data has predicted cell types included
https://tinyurl.com/4yp3dpuc

# read in the data in with LoadSeuratRds
seurat_combined = LoadSeuratRds("study3_seurat_combined.Rds") # 13939 cells across 13 samples, two main groups (CH=clonal hematopoiesis vs contols)



# ------ study data 4
Study title: Single-cell RNA sequencing of peripheral blood reveals immune cell signatures in Alzheimers disease
Paper: https://pubmed.ncbi.nlm.nih.gov/34447367/
original GEO source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181279

# you can download the Rds data here - this data has predicted cell types included
https://tinyurl.com/3m7s3fyk

# read in the data in with LoadSeuratRds
seurat_combined = LoadSeuratRds("study4_seurat_combined.Rds") # 36849 cells across 5 samples, two main groups (Alzheimers vs control)
