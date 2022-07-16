# Install Packages
library(Seurat)
library(dplyr)
library(patchwork)

oncolytic <- Read10X("data/fuson")

# Create Seurat Object
oncolytic_seu <- CreateSeuratObject(oncolytic, min.cells = 3)

#Filtering
VlnPlot(oncolytic_seu, c('nCount_RNA', 'nFeature_RNA'))
FeatureScatter(oncolytic_seu, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

# Add number of genes per UMI for each cell to metadata (This metric with give us an idea of the complexity of our dataset)
oncolytic_seu$log10GenesPerUMI <- log10(oncolytic_seu$nFeature_RNA) / log10(oncolytic_seu$nCount_RNA)

# mtRNA Filtering (This MT pattern works on Homo sapiens datasets)
rownames_oncolytic <- rownames(oncolytic_seu)
rownames_oncolytic[grepl('MT', rownames_oncolytic)]

# Compute percent mt ratio
oncolytic_seu$mtRatio <- PercentageFeatureSet(object = oncolytic_seu, pattern = '^MT-')
oncolytic_seu$mtRatio <- oncolytic_seu@meta.data$mtRatio / 100

# 
# # Create metadata dataframe
# metadata <- oncolytic_seu@meta.data
# 
# # Add cell IDs to metadata
# metadata$cells <- rownames(metadata)
# 
# # Rename columns
# metadata <- metadata %>%
#         dplyr::rename(seq_folder = orig.ident,
#                       nUMI = nCount_RNA,
#                       nGene = nFeature_RNA)
# 
# # Create sample column
# install.packages("stringr")
# library(stringr)
# metadata$sample <- NA
# metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
# metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"
# 
# # Add metadata back to Seurat object
# oncolytic_seu@meta.data <- metadata
# 
# # Create .RData object to load at any time
# save(oncolytic_seu, file="data/merged_filtered_seurat.RData")
# 
# # Visualize the number of cell counts per sample
# install.packages("ggplot2")
# library(ggplot2)
# metadata %>% 
#         ggplot2(aes(x=sample, fill=sample)) + 
#         geom_bar() +
#         theme_classic() +
#         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#         theme(plot.title = element_text(hjust=0.5, face="bold")) +
#         ggtitle("NCells")

VlnPlot(oncolytic_seu, c("nCount_RNA", "nFeature_RNA", "mtRatio"))

oncolytic_seu@meta.data$mtRatio <- NULL

oncolytic_seu <- subset(oncolytic_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 3500)

#Normalize
oncolytic_seu <- NormalizeData(oncolytic_seu)
range(oncolytic_seu@assays$RNA@counts)
range(oncolytic_seu@assays$RNA@data)

#Load and Feature Selection
oncolytic_seu <- FindVariableFeatures(oncolytic_seu, nfeatures = 2000)
top20 <- head(VariableFeatures(oncolytic_seu), 20)
plot_features <- VariableFeaturePlot(oncolytic_seu)
plot_features_label <- LabelPoints(plot = plot_features, points = top20, repel = T)
plot_features_label

# Scale Data
oncolytic_seu <- ScaleData(oncolytic_seu, verbose = FALSE)

#PCA
oncolytic_seu <- RunPCA(oncolytic_seu, verbose = FALSE, npcs = 30)
ElbowPlot(oncolytic_seu, ndims = 30)
DimPlot(oncolytic_seu, reduction = "pca")
DimHeatmap(oncolytic_seu, dims = 1:6, balanced = T)

#UMAP and t-SNE
oncolytic_seu <- RunTSNE(oncolytic_seu, verbose =F, dims = 1:10)
DimPlot(oncolytic_seu, reduction = "tsne")
oncolytic_seu <- RunUMAP(oncolytic_seu, verbose = FALSE, dims = 1:10)
DimPlot(oncolytic_seu, reduction = "umap")

#Clustering
oncolytic_seu <- FindNeighbors(oncolytic_seu, dims = 1:10)
oncolytic_seu <- FindClusters(oncolytic_seu, resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0))
g1 <- DimPlot(oncolytic_seu, reduction = "umap", label = T, group.by = "RNA_snn_res.0.2")
g2 <- DimPlot(oncolytic_seu, reduction = "umap", label = T, group.by = "RNA_snn_res.0.4")
g3 <- DimPlot(oncolytic_seu, reduction = "umap", label = T, group.by = "RNA_snn_res.0.5")
g4 <- DimPlot(oncolytic_seu, reduction = "umap", label = T, group.by = "RNA_snn_res.0.6")
g5 <- DimPlot(oncolytic_seu, reduction = "umap", label = T, group.by = "RNA_snn_res.0.8")
g6 <- DimPlot(oncolytic_seu, reduction = "umap", label = T, group.by = "RNA_snn_res.1")
g1 + g2 + g3 + g4 + g5 + g6

#FindMarkers
Idents(oncolytic_seu) <- "RNA_snn_res.0.8"
marker_table_oncolytic_seu <- FindAllMarkers(oncolytic_seu)

# Visualize
DimPlot(oncolytic_seu, reduction = "umap", label = T)

# Top 20 markers
marker_table_oncolytic_seu %>%
        group_by(cluster) %>%
        top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(oncolytic_seu, features = top20$gene) + NoLegend()

#Visualize Genes
FeaturePlot(oncolytic_seu, reduction = "umap", features = c("Hist1h2ap", "Cxcl14", "C1qa", "Apoe", "Igfbp7", "Hist1h1b", "S100a9", "Il1rl1", "Nkg7", "Cd74", "Stmn1", "Ero1l", "Plac8", "Lxn", "Mgp", "Dcn"))
VlnPlot(oncolytic_seu, features = c("Hist1h2ap", "Cxcl14", "C1qa", "Apoe", "Igfbp7", "Hist1h1b", "S100a9", "Il1rl1", "Nkg7", "Cd74", "Stmn1", "Ero1l", "Plac8", "Lxn", "Mgp", "Dcn"))

#Annotation


