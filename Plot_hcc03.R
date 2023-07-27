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
library(rio)


work_dir <- "./PATH_TO_WORKING_DIR"
output_dir <- './PATH_TO_OUTPUT_DIR'
raw_data_path <- "./PATH_TO_HCC-3_RAW"
ST_result_path <- "./PATH_TO_STLEARN_RESULT"
SCENIC_output_dir <- './PATH_TO_SCENIC_OUTPUT/'
setwd(work_dir)

st_hcc03 <- Load10X_Spatial(raw_data_path)
data_SME_HCC03_identity <- read.csv(sprintf("%s/data_SME_hcc03_identity.csv", ST_result_path), header = T)

######################## Pre-processing #########################
# Assign cell population info generated with StLearn
new_meta <- st_hcc03@meta.data
new_meta$new.ident <- data_SME_HCC03_identity$X_pca_kmeans[which(rownames(new_meta) %in% data_SME_HCC03_identity$X)]
st_hcc03@meta.data <- new_meta

## Visualize clustering
Idents(st_hcc03) <- 'new.ident'
SpatialDimPlot(st_hcc03, cells.highlight = CellsByIdentities(object = st_hcc03, idents = c(0, 1, 2, 3, 4, 5)), 
               facet.highlight = T, ncol = 3)

## QC and filtering
st_hcc03[['percent_mt']] <- PercentageFeatureSet(st_hcc03, pattern = '^MT-')
st_hcc03 <- subset(st_hcc03, subset = nFeature_Spatial > 100 & percent_mt <  40)
# removing spots with zero counts
st_hcc03 <- subset(st_hcc03, subset = nCount_Spatial > 0)
# normalizing with SCTransform
st_hcc03 <- SCTransform(st_hcc03, assay = "Spatial", vars.to.regress = 'percent_mt',verbose = F)

# FindSpatiallyVariables search for features without pre-annotation
# "markvariogram" finds genes whose expression depends on spatial location
st_hcc03 <- FindVariableFeatures(st_hcc03, selection.method = 'vst', nfeatures = 2000)
#top10_varfeat <- head(VariableFeatures(st_hcc03), 10)
# "IGKC" "IGLC2" "AL627171.2" "IGLC3" "IGHA1" "IGHG1" "JCHAIN" "IGHA2" "IGHG4" "IGLC1"  

st_hcc03 <- FindSpatiallyVariableFeatures(st_hcc03, assay = "SCT", features = VariableFeatures(st_hcc03)[1:500], 
                                          selection.method = "markvariogram")
#top_features <- head(SpatiallyVariableFeatures(st_hcc03, selection.method = "markvariogram"), 10)
# "FTL" "AL627171.2" "FTH1" "APOA1" "MT-CO3" "MYL9" "GPX3" "CTSK" "FABP4" "HP"   

##### SCALING THE DATA
all_genes <- rownames(st_hcc03)
st_hcc03 <- ScaleData(st_hcc03, features = all_genes)

##### DIMENSIONALITY REDUCTION, CLUSTERING, VISUALIZATION
st_hcc03 <- RunPCA(st_hcc03, assay = "SCT", verbose = F)
print(st_hcc03[['pca']])

# run UMAP for dimensional reduction
st_hcc03 <- RunUMAP(st_hcc03, reduction = "pca", dims = 1:50)
DimPlot(st_hcc03, reduction = "umap", label = T, split.by = 'new.ident')

# Subset out cluster5 (little spots)
st_hcc03 <- subset(st_hcc03, idents = c(0, 2, 3))
SpatialDimPlot(st_hcc03, cells.highlight = CellsByIdentities(object = st_hcc03, idents = c(0, 2, 3)), 
               facet.highlight = T, ncol = 3)

# Find highly expressed genes in each cluster
hcc03_markers <- FindAllMarkers(st_hcc03, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
markers_df <- as.data.frame(hcc03_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 10, order_by = avg_log2FC))
write.csv(markers_df,
          sprintf("%s/hcc03_markers_stlearn.csv",output_dir),
          row.names = FALSE)

#0 MT1X CYP2E1 ORM1 CYP3A4 ORM2 - tumor + immune
#1 AL627171.2 ALB MT-ND5 MALAT1 NEAT1 - tumor
#4 FABP4 FTH1 CTSK CRYAB AKR1B1 - FAt cells (adipocyte)
#2 IGKC IGHG4 IGLC2 MGP TAGLN - B cells
#3 IGLC2 IGKC IGHA1 IGHA2 IGLC3 - B cells

SpatialFeaturePlot(object = st_hcc03, features = c('CYP1A2','MALAT1','MARCO','IGLC2'), 
                   alpha = c(0.1, 1), ncol = 2)

### Annotate each cell cluster
current.cluster.ids <- c(0,2,3)
new.cluster.ids <- c("Tumor (MT1X)",  "B cell-1 (IGKC)", "B cell-2 (IGLC2)")
st_hcc03@meta.data$seurat_clusters <- plyr::mapvalues(x = st_hcc03@meta.data$new.ident, from = current.cluster.ids, to = new.cluster.ids)

SpatialDimPlot(st_hcc03, group.by = "seurat_clusters") + 
  theme(text = element_text(size = 20)) &
  guides(fill=guide_legend(override.aes=list(size=5)))

DoHeatmap(object = st_hcc03, 
          features = markers_df$gene,
          group.by = 'seurat_clusters')  + 
  theme(text = element_text(size = 20))


write.csv(st_hcc03@meta.data,
          sprintf("%s/st_hcc03_metadata.csv",output_dir),
          row.names = TRUE)


annotate_meta$annotate.ident <- combined_meta$seurat_clusters[which(rownames(annotate_meta) %in% rownames(combined_meta))]
st_hcc03@meta.data <- annotate_meta

SpatialDimPlot(st_hcc03, group.by = 'annotate.ident')
SpatialDimPlot(st_hcc03, group.by = 'seurat_clusters')

save(st_hcc03, 
     file = sprintf("%s/HCC03_UMAPreduced.Rda",output_dir))
#############################  TF module loading  #############################
auc_mtx = (read.table(sprintf('./%s/auc_mtx.csv', SCENIC_output_dir), header = TRUE, row.names = 1, 
                      stringsAsFactors = FALSE, sep = ',')) #load pyscenic output

auc_sub <- auc_mtx[c('XBP1_', 'JUN_', 'FOS_', 'EGR1_',
                     'FLI1_', 'ETS1_', 'NR2F2_'), ]

auc_sub <- t(auc_sub)
colnames(auc_sub) <- c('XBP1_Activity', 'JUN_Activity', 'FOS_Activity', 'EGR1_Activity', 
                       'FLI1_Activity', 'ETS1_Activity', 'NR2F2_Activity')


st_hcc03@meta.data <- merge(st_hcc03@meta.data, auc_sub,
                            by = 'row.names', all = TRUE)

rownames(st_hcc03@meta.data) <- st_hcc03@meta.data$Row.names

#-- violin
CAF_list <- list(c('COL6A3', 'COL3A1', 'COL1A1', 'COL1A2' ,'VIM', 'DCN'))
st_hcc03 <- AddModuleScore(object = st_hcc03, features = CAF_list, name = "CAF_score")
SpatialFeaturePlot(object = st_hcc03, features = "CAF_score1", alpha = c(0.1, 1))

BCellList <- list(c('CD19', 'CD79A', 'CD40', 'CR2', 'FCER2A',
                    'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                    'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'))
st_hcc03 <- AddModuleScore(object = st_hcc03, features = BCellList, name = "Bcell_score")
SpatialFeaturePlot(object = st_hcc03, features = "Bcell_score1", alpha = c(0.1, 1))

plasemaCellList <- list(c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                          'IGHM', 'IGHD', 'IGHG3','IGHG1',
                          'IGKC','IGLC1','IGLC2', 'IGHA1'))
st_hcc03 <- AddModuleScore(object = st_hcc03, features = plasemaCellList, name = "Plasma_score")
SpatialFeaturePlot(object = st_hcc03, features = "Plasma_score1", alpha = c(0.1, 1))

angioList <- list(c('PECAM1', 'VWF', 'ENG', 'CDH5', 'CD34'))
st_hcc03  <- AddModuleScore(object = st_hcc03 , features = angioList, name = "VasScore")
SpatialFeaturePlot(object = st_hcc03 , features = "VasScore1", alpha = c(0.1, 1))


st_hcc03_sub <- subset(x = st_hcc03,  idents = c(2, 3, 4))

SpatialDimPlot(st_hcc03_sub, group.by = "seurat_clusters")

################################     TF module         ################################## 
#-- spatial
SpatialFeaturePlot(object = st_hcc03_sub, 
                   features = c('XBP1_Activity','ETS1_Activity',
                                'JUN_Activity', 'FOS_Activity' ) , 
                   alpha = c(0.1, 1), ncol=2) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


VlnPlot(object = st_hcc03_sub, features = c('XBP1_Activity', 'ETS1_Activity',
                                            'FOS_Activity',  'JUN_Activity'),
        split.by = 'seurat_clusters', ncol=2) & 
  stat_compare_means(size=7)  & 
  theme(text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')


SpatialFeaturePlot(object = st_hcc03_sub, 
                   features = c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                                'IGHM', 'IGHG3','IGHG1',
                                'IGKC','IGLC1','IGLC2', 'IGHA1'), 
                   alpha = c(0.2, 1), ncol=6) &  theme(legend.position='right') 


########################## Angiokine ###########################

angiokine <- c('ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
               'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
               'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
               'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')

angiokine_hcc03 <- c( 'ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
                    'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
                    'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
                    'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')


VlnPlot(object = st_hcc03, features = angiokine,
        split.by = 'seurat_clusters', ncol=6, pt.size=0) & 
  stat_compare_means(size=4)  & 
  theme(text = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')

################################          Spatial           ########################################
SpatialFeaturePlot(object = st_hcc03, 
                   features = c('PRDM1', 'CXCR4', 'COL6A3', 'COL3A1', 'COL1A1',
                                'CD27', 'SDC1','COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=5) & 
  theme(legend.position='right')


SpatialFeaturePlot(object = st_hcc03, 
                   features = c('PRDM1', 'CXCR4', 'CD27', 'SDC1') , 
                   alpha = c(0.2, 1), ncol=2) & 
  theme(legend.position='right') & scale_fill_continuous(type = "viridis")


SpatialFeaturePlot(object = st_hcc03, 
                   features = c('COL6A3', 'COL3A1', 'COL1A1',
                                'COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=3) & 
  theme(legend.position='right') & scale_fill_continuous(type = "viridis")


SpatialFeaturePlot(object = st_hcc03, 
                   features = c('EPCAM', 'THY1', 'KRT19', 'PROM1',
                                'ALDH1A1', 'CD24', 'ANPEP','CD44',
                                'ICAM1', 'CD47', 'SOX9', 'ALDH2', 'ABCG2'), 
                   alpha = c(0.1, 1), ncol=4) & scale_fill_continuous(low = "white", high = "blue")

SpatialFeaturePlot(object = st_hcc03, 
                   features = c('CD19', 'CD79A', 'CD40', 'CR2', 'FCER2A',
                                'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                                'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'), 
                   alpha = c(0.1, 1), ncol=4) & scale_fill_continuous(low = "white", high = "blue")


SpatialFeaturePlot(object = st_hcc03, 
                   features = c('COL6A3', 'COL3A1', 'COL1A1',
                                'COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=3) & 
  theme(legend.position='right')

SpatialFeaturePlot(object = st_hcc03, 
                   features = c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                                'IGHM', 'IGHG3','IGHG1',
                                'IGKC','IGLC1','IGLC2', 'IGHA1'), 
                   alpha = c(0.2, 1), ncol=4) 