library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(EnhancedVolcano)
library(ggpubr)
library(stringr)
library(ggplot2)
library(scales)

work_dir <- "./PATH_TO_WORKING_DIR"
output_dir <- './PATH_TO_OUTPUT_DIR'
raw_data_path <- "./PATH_TO_HCC-2_RAW"
ST_result_path <- "./PATH_TO_STLEARN_RESULT"
SCENIC_output_dir <- './PATH_TO_SCENIC_OUTPUT/'
setwd(work_dir)

st_hcc02 <- Load10X_Spatial(raw_data_path)
data_SME_HCC02_identity <- read.csv(sprintf("%s/data_SME_hcc02_identity.csv", ST_result_path), header = T)

######################## Pre-processing #########################
# Assign cell population info generated with StLearn
new_meta <- st_hcc02@meta.data
new_meta$new.ident <- data_SME_HCC02_identity$X_pca_kmeans[which(rownames(new_meta) %in% data_SME_HCC02_identity$X)]

st_hcc02@meta.data <- new_meta

## Visualize clustering
Idents(st_hcc02) <- 'new.ident'
SpatialDimPlot(st_hcc02, cells.highlight = CellsByIdentities(object = st_hcc02, idents = c(0, 1, 2, 3, 4, 5)), 
               facet.highlight = T, ncol = 3)

## QC and filtering
st_hcc02[['percent_mt']] <- PercentageFeatureSet(st_hcc02, pattern = '^MT-')
st_hcc02 <- subset(st_hcc02, subset = nFeature_Spatial > 100 & percent_mt < 40)
# removing spots with zero counts
st_hcc02 <- subset(st_hcc02, subset = nCount_Spatial > 0)
# normalizing with SCTransform
st_hcc02 <- SCTransform(st_hcc02, assay = "Spatial", vars.to.regress = 'percent_mt',verbose = F)

# FindSpatiallyVariables search for features without pre-annotation
# "markvariogram" finds genes whose expression depends on spatial location
st_hcc02 <- FindVariableFeatures(st_hcc02, selection.method = 'vst', nfeatures = 2000)
top10_varfeat <- head(VariableFeatures(st_hcc02), 10)
# "SPP1"   "IGLC1" "IGLC3"  "MMP7"   "IGHG4" "IGHM"   "AL627171.2" "MMP9"   "IGLC7" "AC240274.1"

st_hcc02 <- FindSpatiallyVariableFeatures(st_hcc02, assay = "SCT", features = VariableFeatures(st_hcc02)[1:500], 
                                          selection.method = "markvariogram")
top_features <- head(SpatiallyVariableFeatures(st_hcc02, selection.method = "markvariogram"), 10)
# "MALAT1" "ALB"    "SPP1"   "HP"     "APOA1"  "FGA"    "SAA1"   "AMBP"   "ORM1"   "SAA2"  

##### SCALING THE DATA
all_genes <- rownames(st_hcc02)
st_hcc02 <- ScaleData(st_hcc02, features = all_genes)

##### DIMENSIONALITY REDUCTION, CLUSTERING, VISUALIZATION
st_hcc02 <- RunPCA(st_hcc02, assay = "SCT", verbose = F)
print(st_hcc02[['pca']])

# run UMAP for dimensional reduction
st_hcc02 <- RunUMAP(st_hcc02, reduction = "pca", dims = 1:50)
DimPlot(st_hcc02, reduction = "umap", label = T)

st_hcc02 <- subset(st_hcc02, idents = c(0, 1, 2, 4))
SpatialDimPlot(st_hcc02, cells.highlight = CellsByIdentities(object = st_hcc02, idents = c(0, 1, 2, 4)), 
               facet.highlight = T, ncol = 3)

# Find highly expressed genes in each cluster
hcc02_markers <- FindAllMarkers(st_hcc02, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
markers_df <- as.data.frame(hcc02_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 10, order_by = avg_log2FC))
write.csv(markers_df,
          "Long_result/clusters_marker_genes_stlearn/hcc02_markers_stlearn.csv", 
          row.names = FALSE)

SpatialFeaturePlot(object = st_hcc02, features = c('ALB','LRRC75A','PLA2G1B','IGKC','COL1A1'), 
                   alpha = c(0.1, 1), ncol = 3)

# 1: 
# 2: 
### Annotate each cell cluster
current.cluster.ids <- c(0, 1, 2, 4)
new.cluster.ids <- c("Myeloid cell", "Immune", "Tumor-1 (ALB)", "Tumor-2 (MALAT1)")
st_hcc02@meta.data$seurat_clusters <- plyr::mapvalues(x = st_hcc02@meta.data$new.ident, from = current.cluster.ids, to = new.cluster.ids)

SpatialDimPlot(st_hcc02, group.by = "seurat_clusters") + 
  theme(text = element_text(size = 20)) &
  guides(fill=guide_legend(override.aes=list(size=5)))

DoHeatmap(object = st_hcc02, 
          features = markers_df$gene,
          group.by = 'seurat_clusters')  + 
  theme(text = element_text(size = 20))


write.csv(st_hcc02@meta.data,
          sprintf("%s/st_hcc02_metadata.csv",output_dir),
          row.names = TRUE)



annotate_meta$annotate.ident <- combined_meta$seurat_clusters[which(rownames(annotate_meta) %in% rownames(combined_meta))]
st_hcc02@meta.data <- annotate_meta

SpatialDimPlot(st_hcc02, group.by = 'annotate.ident')

save(st_hcc02, 
     file = sprintf("%s/hcc02_UMAPreduced.Rda",output_dir))


###########################     Manuscript             ##############################
SpatialDimPlot(st_hcc02, group.by = "seurat_clusters")

auc_mtx = (read.table(sprintf('./%s/auc_mtx.csv', SCENIC_output_dir), header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE, sep = ',')) #load pyscenic output

auc_sub <- auc_mtx[c('PAX5_', 'TCF4_', 'STAT1_', 
                     'XBP1_', 'IRF4_', 'FLI1_',
                     'JUN_', 'FOS_', 'EGR1_',
                     'FLI1_', 'ETS1_', 'NR2F2_'),]

auc_sub <- t(auc_sub)
colnames(auc_sub) <- c('PAX5_Activity', 'TCF4_Activity', 'STAT1_Activity',
                       'XBP1_Activity', 'IRF4_Activity', 'FLI1_Activity',
                       'JUN_Activity', 'FOS_Activity', 'EGR1_Activity', 
                       'FLI1_Activity', 'ETS1_Activity', 'NR2F2_Activity')


st_hcc02@meta.data <- merge(st_hcc02@meta.data, auc_sub,
                            by = 'row.names', all = TRUE)

rownames(st_hcc02@meta.data) <- st_hcc02@meta.data$Row.names

#-- violin
CAF_list <- list(c('COL6A3', 'COL3A1', 'COL1A1', 'COL1A2' ,'VIM', 'DCN'))
st_hcc02 <- AddModuleScore(object = st_hcc02, features = CAF_list, name = "CAF_score")
SpatialFeaturePlot(object = st_hcc02, features = "CAF_score1", alpha = c(0.1, 1))

plasemaCellList <- list(c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                          'IGHM', 'IGHD', 'IGHG3','IGHG1',
                          'IGKC','IGLC1','IGLC2', 'IGHA1'))
st_hcc02 <- AddModuleScore(object = st_hcc02, features = plasemaCellList, name = "Plasma_score")
SpatialFeaturePlot(object = st_hcc02, features = "Plasma_score1", alpha = c(0.1, 1))

BCellList <- list(c('CD19', 'CD79A', 'CD40', 'CR2', 'FCER2A',
                    'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                    'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'))
st_hcc02 <- AddModuleScore(object = st_hcc02, features = BCellList, name = "Bcell_score")
SpatialFeaturePlot(object = st_hcc02, features = "Bcell_score1", alpha = c(0.1, 1))

angioList <- list(c('PECAM1', 'VWF', 'ENG', 'CDH5', 'CD34'))
st_hcc02 <- AddModuleScore(object = st_hcc02, features = angioList, name = "VasScore")
SpatialFeaturePlot(object = st_hcc02, features = "VasScore1", alpha = c(0.1, 1))

# subset myeloid and b cell 
st_hcc02_sub <- subset(x = st_hcc02, idents = c(0, 1))


################################     TF module         ################################## 
st_hcc02 <- subset(x = st_hcc02, subset = PAX5_Activity < 0.5)


#-- spatial


SpatialFeaturePlot(object = st_hcc02, 
                   features = c('PAX5_Activity', 'XBP1_Activity', 'JUN_Activity', 'FOS_Activity' ) , 
                   alpha = c(0.1, 1), ncol=2) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


VlnPlot(object = st_hcc02_sub, features = c('XBP1_Activity', 'ETS1_Activity',
                                            'FOS_Activity',  'JUN_Activity'),
        split.by = 'seurat_clusters', ncol=2) & 
  stat_compare_means(size=6)  & 
  theme(text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')


SpatialFeaturePlot(object = st_hcc02, 
                   features = c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                                'IGHM', 'IGHD', 'IGHG3','IGHG1',
                                'IGKC','IGLC1','IGLC2', 'IGHA1'), 
                   alpha = c(0.2, 1), ncol=4) & scale_fill_continuous(low = "white", high = "red") & 
  theme(legend.position='right') 



################################          Spatial           ########################################
SpatialFeaturePlot(object = st_hcc02, 
                   features = c('PRDM1', 'CXCR4', 'COL6A3', 'COL3A1', 'COL1A1',
                                'CD27', 'SDC1','COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=5) & 
  theme(legend.position='right')



SpatialFeaturePlot(object = st_hcc02, 
                   features = c('EPCAM', 'THY1', 'KRT19', 'PROM1',
                                'ALDH1A1', 'CD24', 'ANPEP','CD44',
                                'ICAM1', 'CD47', 'SOX9', 'ALDH2', 'ABCG2'), 
                   alpha = c(0.1, 1), ncol=4) & scale_fill_continuous(low = "white", high = "blue")

SpatialFeaturePlot(object = st_hcc02, 
                   features = c('CD19', 'CD79A', 'CD40', 'CR2', 'FCER2A',
                                'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                                'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'), 
                   alpha = c(0.1, 1), ncol=4) & scale_fill_continuous(low = "white", high = "blue")


SpatialFeaturePlot(object = st_hcc02, 
                   features = c('COL6A3', 'COL3A1', 'COL1A1',
                                'COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=3) & 
  theme(legend.position='right')

SpatialFeaturePlot(object = st_hcc02, 
                   features = c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                                'IGHM', 'IGHD', 'IGHG3','IGHG1',
                                'IGKC','IGLC1','IGLC2', 'IGHA1'), 
                   alpha = c(0.2, 1), ncol=4) 

