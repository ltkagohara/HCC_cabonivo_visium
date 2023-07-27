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
raw_data_path <- "./PATH_TO_HCC-4_RAW"
ST_result_path <- "./PATH_TO_STLEARN_RESULT"
SCENIC_output_dir <- './PATH_TO_SCENIC_OUTPUT/'
setwd(work_dir)

st_hcc04 <- Load10X_Spatial(raw_data_path)
data_SME_HCC04_identity <- read.csv(sprintf("%s/data_SME_hcc04_identity.csv", ST_result_path), header = T)

######################## Pre-processing #########################
# Assign cell population info generated with StLearn

new_meta <- st_hcc04@meta.data
new_meta$new.ident <- data_SME_HCC04_identity$X_pca_kmeans[which(rownames(new_meta) %in% data_SME_HCC04_identity$X)]
st_hcc04@meta.data <- new_meta

## Visualize clustering
Idents(st_hcc04) <- 'new.ident'
SpatialDimPlot(st_hcc04, cells.highlight = CellsByIdentities(object = st_hcc04, idents = c(0, 1, 2, 3, 4, 5, 6, 7,8,9)), 
               facet.highlight = T, ncol = 5)

## QC and filtering
st_hcc04[['percent_mt']] <- PercentageFeatureSet(st_hcc04, pattern = '^MT-')
st_hcc04 <- subset(st_hcc04, subset = nFeature_Spatial > 100 & percent_mt < 40)
# removing spots with zero counts
st_hcc04 <- subset(st_hcc04, subset = nCount_Spatial > 0)
# normalizing with SCTransform
st_hcc04 <- SCTransform(st_hcc04, assay = "Spatial", vars.to.regress = 'percent_mt',verbose = F)

# FindSpatiallyVariables search for features without pre-annotation
# "markvariogram" finds genes whose expression depends on spatial location
st_hcc04 <- FindVariableFeatures(st_hcc04, selection.method = 'vst', nfeatures = 2000)
#top10_varfeat <- head(VariableFeatures(st_hcc04), 10)
# "IGKC"       "IGHM"       "IGLC3"      "COL1A1"     "AL627171.2" "IGLC2"      "CD74"       "MALAT1"     "IGHA1"      "IGHG3"  

st_hcc04 <- FindSpatiallyVariableFeatures(st_hcc04, assay = "SCT", features = VariableFeatures(st_hcc04)[1:500], 
                                          selection.method = "markvariogram")
#top_features <- head(SpatiallyVariableFeatures(st_hcc04, selection.method = "markvariogram"), 10)
# "APOC1"   "CD74"    "ALB"     "MALAT1"  "APOA2"   "APOA1"   "HLA-DRA" "TIMP1"   "FGA"     "COL1A1"     

##### SCALING THE DATA
all_genes <- rownames(st_hcc04)
st_hcc04 <- ScaleData(st_hcc04, features = all_genes)

##### DIMENSIONALITY REDUCTION, CLUSTERING, VISUALIZATION
st_hcc04 <- RunPCA(st_hcc04, assay = "SCT", verbose = F)
print(st_hcc04[['pca']])

# run UMAP for dimensional reduction
st_hcc04 <- RunUMAP(st_hcc04, reduction = "pca", dims = 1:50)
DimPlot(st_hcc04, reduction = "umap", label = T, split.by = 'new.ident')

# Subset out cluster5 (little spots) and cluster3/6/9/4 (dead cells)
st_hcc04 <- subset(st_hcc04, idents = c(0, 1, 2, 7, 8))
SpatialDimPlot(st_hcc04, cells.highlight = CellsByIdentities(object = st_hcc04, idents = c(0, 1, 2, 7, 8)), 
               facet.highlight = T, ncol = 3)

# Find highly expressed genes in each cluster
hcc04_markers <- FindAllMarkers(st_hcc04, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)

markers_df <- as.data.frame(hcc04_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 10, order_by = avg_log2FC))

SpatialDimPlot(st_hcc04, group.by = "seurat_clusters") + 
  theme(text = element_text(size = 20)) &
  guides(fill=guide_legend(override.aes=list(size=5)))

DoHeatmap(object = st_hcc04, 
          features = markers_df$gene,
          group.by = 'seurat_clusters')  + 
  theme(text = element_text(size = 20))

write.csv(markers_df,
          sprintf("%s/hcc04_markers_stlearn.csv",output_dir),
          row.names = FALSE)
#9 YBX1 MT-ATP8 AL627171.2 AL355075.4 MT-ND4L  -dead
#6 AL627171.2 MT-ND5 AL355075.4 MT-ND2 MT-CYB -dead
#3 AL627171.2 AL355075.4 MT-ATP8 MT-ND4L MT-ND5 -dead
#4 MT-ND5 AL627171.2 MT-CYB MT-ND2 -dead

#2 AL031058.1 CST1 APCDD1 NOTUM SLC1A5 - tumor/immune
#1 UGT2B15 ORM1 APOC4 HSD11B1 ITIH2 - tumor/immune
#0 CYP2A6 HAMP FABP1 MALAT1 APOC3 -tumor 
#7 IGKC IGHM CD52 CCL19 CD74 - Immune
#8 COL1A1 COL1A2 COL3A1 AEBP1 CD74 - CAF

### Annotate each cell cluster
current.cluster.ids <- c(0, 1, 2, 7, 8)
new.cluster.ids <- c("Tumor-1 (CYP2A6)", "Tumor-2 (UGT2B15)", "Tumor-3 (CST1)", "Immune", "CAFs")
st_hcc04@meta.data$seurat_clusters <- plyr::mapvalues(x = st_hcc04@meta.data$new.ident, from = current.cluster.ids, to = new.cluster.ids)

SpatialDimPlot(st_hcc04, group.by = "seurat_clusters")


write.csv(st_hcc04@meta.data,
          sprintf("%s/st_hcc04_metadata.csv",output_dir),
          row.names = TRUE)

save(st_hcc04, 
     file = sprintf("%s/hcc01_UMAPreduced.Rda",output_dir))

######################### TF module plotting #######################
SpatialDimPlot(st_hcc04, group.by = "seurat_clusters")

auc_mtx <- (read.table(sprintf('./%s/auc_mtx.csv', SCENIC_output_dir), header = TRUE, row.names = 1, 
                      stringsAsFactors = FALSE, sep = ',')) #load pyscenic output

auc_sub <- auc_mtx[c('PAX5_', 'TCF4_', 'STAT1_', 
                     'IRF4_', 'IRF5_', 'IRF1_',
                     'JUN_', 'FOS_', 'EGR1_',
                     'XBP1_'), ]

auc_sub <- t(auc_sub)
colnames(auc_sub) <- c('PAX5_Activity', 'TCF4_Activity', 'STAT1_Activity',
                       'IRF4_Activity', 'IRF5_Activity', 'IRF1_Activity',
                       'JUN_Activity', 'FOS_Activity', 'EGR1_Activity',
                       'XBP1_Activity')


st_hcc04@meta.data <- merge(st_hcc04@meta.data, auc_sub,
                            by = 'row.names', all = TRUE)

rownames(st_hcc04@meta.data) <- st_hcc04@meta.data$Row.names

#-- violin
CAF_list <- list(c('COL6A3', 'COL3A1', 'COL1A1', 'COL1A2' ,'VIM', 'DCN'))
st_hcc04 <- AddModuleScore(object = st_hcc04, features = CAF_list, name = "CAF_score")
SpatialFeaturePlot(object = st_hcc04, features = "CAF_score1", alpha = c(0.1, 1))

BCellList <- list(c('CD19', 'CD79A', 'CD40', 'CR2', 'FCER2A',
                    'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                    'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'))
st_hcc04 <- AddModuleScore(object = st_hcc04, features = BCellList, name = "Bcell_score")
SpatialFeaturePlot(object = st_hcc04, features = "Bcell_score1", alpha = c(0.1, 1))


plasemaCellList <- list(c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                          'IGHM', 'IGHD', 'IGHG3','IGHG1',
                          'IGKC','IGLC1','IGLC2', 'IGHA1'))
st_hcc04 <- AddModuleScore(object = st_hcc04, features = plasemaCellList, name = "Plasma_score")
SpatialFeaturePlot(object = st_hcc04, features = "Plasma_score1", alpha = c(0.1, 1))

angioList <- list(c('PECAM1', 'VWF', 'ENG', 'CDH5', 'CD34'))
st_hcc04 <- AddModuleScore(object = st_hcc04, features = angioList, name = "VasScore")
SpatialFeaturePlot(object = st_hcc04, features = "VasScore1", alpha = c(0.1, 1))

# 7: immune    8: CAFs
st_hcc04_sub <- subset(x = st_hcc04,  idents = c(7, 8))
st_hcc04_immune<- subset(x = st_hcc04,  idents = 7)


################################     TF module         ################################## 
#-- spatial
SpatialFeaturePlot(object = st_hcc04, 
                   features = c('PAX5_Activity', 'IRF1_Activity',
                                'FOS_Activity',  'JUN_Activity') , 
                   alpha = c(0.1, 1), ncol=2) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


VlnPlot(object = st_hcc04_sub, features = c('PAX5_Activity', 'IRF1_Activity',
                                            'FOS_Activity',  'JUN_Activity'),
        split.by = 'seurat_clusters', ncol=2) & 
  stat_compare_means(size=6)  & 
  theme(text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')


#Figure HCC 14 B
SpatialFeaturePlot(object = st_hcc04_sub, 
                   features = c('PAX5_Activity', 'IRF1_Activity',
                                'FOS_Activity',  'JUN_Activity') , 
                   alpha = c(0.1, 1), ncol=2) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


SpatialFeaturePlot(object = st_hcc04, 
                   features = c('XBP1_Activity') , 
                   alpha = c(0.1, 1), ncol=1) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")



########################## Angiokine ###########################

angiokine <- c('ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
               'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
               'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
               'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')

angiokine_hcc04 <- c( 'ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
                    'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
                    'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
                    'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')


VlnPlot(object = st_hcc04, features = angiokine,
        split.by = 'seurat_clusters', ncol=6, pt.size=0) & 
  stat_compare_means(size=4)  & 
  theme(text = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')


################################          Spatial           ########################################


SpatialFeaturePlot(object = st_hcc04, 
                   features = c('CD79A', 'CD19', 'COL6A3', 'COL3A1', 'COL1A1',
                                'CD22', 'CIITA','COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=5) & 
  theme(legend.position='right')


SpatialFeaturePlot(object = st_hcc04, 
                   features = c('CD79A', 'CD19','CD22', 'CIITA') , 
                   alpha = c(0.2, 1), ncol=2) & 
  theme(legend.position='right') 

SpatialFeaturePlot(object = st_hcc04, 
                   features = c('COL6A3', 'COL3A1', 'COL1A1',
                                'COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=3) & 
  theme(legend.position='right') 


SpatialFeaturePlot(object = tumor_region, 
                   features = c('IGHG3', 'IGHG1', 'IGLC1', 'IGLC2',
                                'CCL21', 'CCL19') , 
                   alpha = c(0.1, 1), ncol=4) & scale_fill_continuous(low = "white", high = "blue")











