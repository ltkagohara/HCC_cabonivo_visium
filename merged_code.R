library(Seurat)
library(dplyr)
library(cowplot)
library(harmony)
library(devtools)
library(SummarizedExperiment)
library(Matrix.utils)
library(DESeq2)
library(RColorBrewer)
library(EnhancedVolcano)

work_dir <- "./PATH_TO_WORKING_DIR"
setwd(work_dir)
output_dir <- './PATH_TO_OUTPUT_DIR'

HCC1_path <- './PATH_TO_HCC1_DATA'
load(HCC1_path)

HCC2_path <- './PATH_TO_HCC2_DATA'
load(HCC2_path)

HCC3_path <- './PATH_TO_HCC3_DATA'
load(HCC3_path)

HCC4_path <- './PATH_TO_HCC4_DATA'
load(HCC4_path)

HCC5_path <- './PATH_TO_HCC5_DATA'
load(HCC5_path)

HCC6_path <- './PATH_TO_HCC6_DATA'
load(HCC6_path)

HCC7_path <- './PATH_TO_HCC7_DATA'
load(HCC7_path)

SpatialDimPlot(st_hcc01, group.by = "seurat_clusters")
st_hcc02$orig.ident <- 'HCC_1'

SpatialDimPlot(st_hcc02, group.by = "seurat_clusters")
st_hcc08$orig.ident <- 'HCC_2'

SpatialDimPlot(st_hcc03, group.by = "seurat_clusters")
st_hcc13$orig.ident <- 'HCC_3'

SpatialDimPlot(st_hcc04, group.by = "seurat_clusters")
st_hcc14$orig.ident <- 'HCC_4'

SpatialDimPlot(st_hcc05, group.by = "seurat_clusters")
st_hcc04$orig.ident <- 'HCC_5'

SpatialDimPlot(st_hcc06, group.by = "seurat_clusters")
st_hcc09$orig.ident <- 'HCC_6'

SpatialDimPlot(st_hcc07, group.by = "seurat_clusters")
st_hcc20$orig.ident <- 'HCC_7'

data_merged <- merge(st_hcc01, 
                     y = c(st_hcc02,st_hcc03,st_hcc04,st_hcc05,st_hcc06, st_hcc07), 
                     add.cell.ids = c("HCC01","HCC02","HCC03","HCC04","HCC05","HCC06", "HCC07"))

data_merged <- subset(data_merged, subset = nCount_Spatial > 0)

#table(tumor.all$orig.ident)
#Assay data with 36601 features for 10318 cells
data.list <- SplitObject(object = data_merged, split.by = "orig.ident")
for (i in 1:length(data.list)) {
  data.list[[i]] <- SCTransform(data.list[[i]], verbose = FALSE, assay = 'Spatial')
}
tumor.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data_merged <- merge(data.list[[1]], 
                     y = data.list[2:length(data.list)],  
                     project = "HCC_All", 
                     merge.data = TRUE)
VariableFeatures(data_merged) <- tumor.features

new.meta <- data_merged@meta.data
new.meta$response <- 'Responder'
new.meta$response[new.meta$orig.ident == c("HCC_5")] <- "Non-responder"
new.meta$response[new.meta$orig.ident == c("HCC_6")] <- "Non-responder"
new.meta$response[new.meta$orig.ident == c("HCC_7")] <- "Non-responder"
data_merged@meta.data <- new.meta


data_merged <- RunPCA(object = data_merged, assay = "SCT", npcs = 50)

# t-SNE and Clustering
data_merged <- RunUMAP(data_merged, reduction = "pca", dims = 1:20)
data_merged <- FindNeighbors(data_merged, reduction = "pca", dims = 1:20)
data_merged <- FindClusters(data_merged, resolution = 0.5)

data_merged <- RunHarmony(data_merged, group.by.vars = "orig.ident", assay.use = "SCT", plot_convergence = T)
harmony_embeddings <- Embeddings(data_merged, 'harmony')

p1 <- DimPlot(object = data_merged, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = data_merged, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)

data_merged <- data_merged %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DefaultAssay(data_merged) <- "SCT"

p3 <- DimPlot(object = data_merged, reduction = "umap", pt.size = .1, group.by = "orig.ident")
p4 <- DimPlot(object = data_merged, reduction = "umap", pt.size = .1, group.by = "seurat_clusters", label = T)
plot_grid(p1,p2)

DimPlot(object = data_merged, reduction = "umap", pt.size = .1, group.by = "response")
DimPlot(object = data_merged, reduction = "umap", pt.size = .1, group.by = "orig.ident")

############################  Antigen Presentation ######################## (add all module scores here)
# Antigen presentation pathway
agpres_genes <- list(c('B2M','CALR','CANX','CD4','CD74','CD8A','CD8B','CIITA','CREB1','CTSB','CTSL','CTSS','HLA-A','HLA-B','HLA-C','HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQA2','HLA-DQB1','HLA-DRA','HLA-DRB1','HLA-DRB3','HLA-DRB4','HLA-DRB5','HLA-E','HLA-F','HLA-G','HSP90AA1','HSP90AB1','HSPA1A','HSPA1B','HSPA1L','HSPA2',
                       'HSPA4','HSPA5','HSPA6','HSPA8','IFI30','IFNA1','IFNA10','IFNA13','IFNA14','IFNA16',
                       'IFNA17','IFNA2','IFNA21','IFNA4','IFNA5','IFNA6','IFNA7','IFNA8','KIR2DL1','KIR2DL2',
                       'KIR2DL3','KIR2DL4','KIR2DL5A','KIR2DS1','KIR2DS3','KIR2DS4','KIR2DS5','KIR3DL1','KIR3DL2','KIR3DL3','KLRC1','KLRC2','KLRC3','KLRC4','KLRD1','LGMN','LTA','NFYA','NFYB','NFYC','PDIA3','PSME1','PSME2','PSME3','RFX5','RFXANK','RFXAP','TAP1','TAP2','TAPBP'))

data_merged <- AddModuleScore(object = data_merged, features = agpres_genes, name = 'AgPres')

#CD39->ENPD1, aSMA->ACTA2
panCAFList <- list(c('LUM', 'DCN', 'COL1A1', 'VIM', 'ENTPD1', 
                     'ACTA2', 'PDPN', 'COL1A2', 'SERPINF2'))
data_merged <- AddModuleScore(object = data_merged, features = panCAFList, name = "panCAFScore")

iCAFList <- list(c('CXCL1', 'CXCL2', 'CCL2', 'CXCL12', 'PDGFRA', 
                   'CFD', 'LMNA', 'DPT', 'HAS1', 'HAS2'))
data_merged <- AddModuleScore(object = data_merged, features = iCAFList, name = "iCAFScore")

myCAFList <- list(c('ACTA2', 'TAGLN', 'MYL9', 'TPM2', 'MMP11', 
                    'POSTN', 'HOPX', 'TWIST1', 'SOX4', 'ZUB1'))
data_merged <- AddModuleScore(object = data_merged, features = myCAFList, name = "myCAFScore")


data_merged <- PrepSCTFindMarkers(data_merged, assay = "SCT", verbose = TRUE)
data_merged_markers <- FindAllMarkers(data_merged,  grouping.var = "orig.ident")
markers_df <- as.data.frame(data_merged_markers %>% group_by(cluster) %>% slice_max(n=20, order_by = avg_log2FC))


DefaultAssay(data_merged) <- "SCT"


DoHeatmap(object = data_merged, 
          features = markers_df$gene,  #cluster5 
          angle = 45, hjust = 0)


current.cluster.ids <- c(0, 1, 2, 3, 4, 5,
                         6, 7, 8, 9, 10, 11,
                         12, 13, 14, 15)

#4,5,8,10
DoHeatmap(object = data_merged, 
          features = markers_df$gene,  #cluster5 
          group.by = 'seurat_clusters', size=4, angle = 45, hjust = 0, label = FALSE,)

DimPlot(object = data_merged, reduction = "umap", pt.size = .1, group.by = "seurat_clusters")

new.cluster.ids <- c("Tumor (PRSS1)", "Tumor (SAA1)", "Tumor (MALAT1)", "Tumor (CYP2E1)", "CAFs", "Immune (IGLC2)", 
                     "Tumor (APOA1)", "Tumor (MT-ND2)" , "Immune (IGHG3)", "Tumor (AFP)", "Macrophage", "Tumor (HP)",
                     "Tumor (HBA2)", "Tumor (MMP7)", "Tumor (H19)", "Tumor (TF)")
data_merged@meta.data$seurat_clusters <- plyr::mapvalues(x = data_merged@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

write.csv(data_merged@meta.data,
          sprintf("%s/merged_metadata.csv",output_dir), 
          row.names = TRUE)


DoHeatmap(object = data_merged, 
          features = markers_df$gene,  #cluster5 
          group.by = 'seurat_clusters', size=4, angle = 45, hjust = 0, label = FALSE,)



# creating a palette of colors that fits 15 clusters
colorcount <- length(unique(data_merged$SCT_snn_res.0.5))
myPalette <- colorRampPalette(brewer.pal(colorcount, "Paired"))
colorValue <- myPalette(colorcount)


colorMap <- c("Tumor (PRSS1)" = colorValue[1], "Tumor (SAA1)" = colorValue[2], "Tumor (MALAT1)" = colorValue[3], "Tumor (CYP2E1)" = colorValue[4], 
              "CAFs" = colorValue[5], "Immune (IGLC2)" = colorValue[6], "Tumor (APOA1)" = colorValue[7], "Tumor (MT-ND2)" = colorValue[8], 
              "Immune (IGHG3)" = colorValue[9], "Tumor (AFP)" = colorValue[10], "Macrophage" = colorValue[11], "Tumor (HP)" = colorValue[12],
              "Tumor (HBA2)" = colorValue[13], "Tumor (MMP7)" = colorValue[14], "Tumor (H19)" = colorValue[15], "Tumor (TF)" = colorValue[16])


DimPlot(object = data_merged, reduction = "umap", cols = colorMap, pt.size = .1, group.by = "seurat_clusters")

DimPlot(object = data_merged, reduction = "umap", pt.size = .1, group.by = "orig.ident")


###Spatial plot
#HCC1
HCC1_spatial <-SpatialDimPlot(data_merged, group.by = "seurat_clusters", 
                              images = c('slice1'), cols = colorMap) &
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC1_spatial <- HCC1_spatial + labs(title = "HCC1-R") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))


HCC1_spatial

#HCC-02-plot
HCC2_spatial <-  SpatialDimPlot(data_merged, group.by = "seurat_clusters", cols = colorMap,  
                                images = c('slice1_HCC02.1')) &
  guides(fill=guide_legend(override.aes=list(size=5))) 

HCC2_spatial <- HCC2_spatial + labs(title = "HCC2-R") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC2_spatial

#HCC3
HCC3_spatial <-SpatialDimPlot(data_merged, group.by = "seurat_clusters", 
                              images = c('slice1_HCC03.2'), cols = colorMap) &
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC3_spatial <- HCC3_spatial + labs(title = "HCC3-R") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC3_spatial

#HCC4
HCC4_spatial <-SpatialDimPlot(data_merged, group.by = "seurat_clusters", 
                              images = c('slice1_HCC04.3'), cols = colorMap)&
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC4_spatial <- HCC4_spatial + labs(title = "HCC4-R") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC4_spatial

#HCC5
HCC5_spatial <-SpatialDimPlot(data_merged, group.by = "seurat_clusters", 
                              images = c('slice1_HCC05.4'), cols = colorMap)&
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC5_spatial <- HCC5_spatial + labs(title = "HCC5-NR") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC5_spatial

#HCC6
HCC6_spatial <-SpatialDimPlot(data_merged, group.by = "seurat_clusters", 
                              images = c('slice1_HCC06.5'), cols = colorMap)&
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC6_spatial <- HCC6_spatial + labs(title = "HCC6-NR") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC6_spatial

##HCC7
HCC7_spatial <-SpatialDimPlot(data_merged, group.by = "seurat_clusters", 
                              images = c('slice1_HCC07.6'), cols = colorMap) &
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC7_spatial <- HCC7_spatial + labs(title = "HCC7-NR") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))


HCC7_spatial
save(data_merged, 
     file = sprintf("%s/data_combined_UMAPreduced.Rda",output_dir))

############################  Antigen Plot / Tumor plots ###########################
tumor_region <- subset(x = data_merged, idents = c(0, 1, 2, 3, 6, 7, 9,
                                                   11, 12, 13, 14, 15))


VlnPlot(tumor_region, features = c('AgPres1'), group.by = 'response')


SpatialFeaturePlot(tumor_region, features = c('AgPres1'), 
                   images = c('slice1', "slice1_HCC02.1")) 



SpatialFeaturePlot(tumor_region, features = c('AgPres1'), 
                   images = c('slice1_HCC03.2', "slice1_HCC04.3")) 



SpatialFeaturePlot(tumor_region, features = c('AgPres1'), 
                   images = c('slice1_HCC05.4', "slice1_HCC06.5")) 



SpatialFeaturePlot(tumor_region, features = c('AgPres1'), 
                   images = c('slice1_HCC7.6')) 



###Spatial plot for tumor region only
#HCC1
HCC1_spatial <-SpatialDimPlot(tumor_region, group.by = "seurat_clusters", 
                              images = c('slice1'), cols = colorMap) &
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC1_spatial <- HCC1_spatial + labs(title = "HCC1-R") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))


HCC1_spatial

#HCC-02-plot
HCC2_spatial <-  SpatialDimPlot(tumor_region, group.by = "seurat_clusters", cols = colorMap,  
                                images = c('slice1_HCC02.1')) &
  guides(fill=guide_legend(override.aes=list(size=5))) 

HCC2_spatial <- HCC2_spatial + labs(title = "HCC2-R") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC2_spatial

#HCC3
HCC3_spatial <-SpatialDimPlot(tumor_region, group.by = "seurat_clusters", 
                              images = c('slice1_HCC03.2'), cols = colorMap) &
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC3_spatial <- HCC3_spatial + labs(title = "HCC3-R") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC3_spatial

#HCC4
HCC4_spatial <-SpatialDimPlot(tumor_region, group.by = "seurat_clusters", 
                              images = c('slice1_HCC04.3'), cols = colorMap)&
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC4_spatial <- HCC4_spatial + labs(title = "HCC4-R") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC4_spatial

#HCC5
HCC5_spatial <-SpatialDimPlot(tumor_region, group.by = "seurat_clusters", 
                              images = c('slice1_HCC05.4'), cols = colorMap)&
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC5_spatial <- HCC5_spatial + labs(title = "HCC5-NR") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC5_spatial

#HCC6
HCC6_spatial <-SpatialDimPlot(tumor_region, group.by = "seurat_clusters", 
                              images = c('slice1_HCC06.5'), cols = colorMap)&
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC6_spatial <- HCC6_spatial + labs(title = "HCC6-NR") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))

HCC6_spatial

##HCC7
HCC7_spatial <-SpatialDimPlot(tumor_region, group.by = "seurat_clusters", 
                              images = c('slice1_HCC07.6'), cols = colorMap) &
  guides(fill=guide_legend(override.aes=list(size=5)))

HCC7_spatial <- HCC7_spatial + labs(title = "HCC7-NR") + 
  theme(plot.title = element_text(hjust = 0.5, size = 24))


HCC7_spatial

######################## 

