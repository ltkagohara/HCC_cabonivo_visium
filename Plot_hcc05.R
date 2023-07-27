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
raw_data_path <- "./PATH_TO_HCC-5_RAW"
ST_result_path <- "./PATH_TO_STLEARN_RESULT"
SCENIC_output_dir <- './PATH_TO_SCENIC_OUTPUT/'
setwd(work_dir)

st_hcc05 <- Load10X_Spatial(raw_data_path)
data_SME_HCC05_identity <- read.csv(sprintf("%s/data_SME_hcc05_identity.csv", ST_result_path), header = T)
############################################ Pre processing ################################################
new_meta <- st_hcc05@meta.data
new_meta$new.ident <- data_SME_HCC05_identity$X_pca_kmeans[which(rownames(new_meta) %in% data_SME_HCC05_identity$X)]
st_hcc05@meta.data <- new_meta

## Visualize clustering
Idents(st_hcc05) <- 'new.ident'
SpatialDimPlot(st_hcc05, cells.highlight = CellsByIdentities(object = st_hcc05, idents = c(0, 1, 2, 3)), 
               facet.highlight = T, ncol = 2)

## QC and filtering
st_hcc05[['percent_mt']] <- PercentageFeatureSet(st_hcc05, pattern = '^MT-')
st_hcc05 <- subset(st_hcc05, subset = nFeature_Spatial > 100 & percent_mt < 40)
# removing spots with zero counts
st_hcc05 <- subset(st_hcc05, subset = nCount_Spatial > 0)
# normalizing with SCTransform
st_hcc05 <- SCTransform(st_hcc05, assay = "Spatial", vars.to.regress = 'percent_mt',verbose = F)

# FindSpatiallyVariables search for features without pre-annotation
# "markvariogram" finds genes whose expression depends on spatial location
st_hcc05 <- FindVariableFeatures(st_hcc05, selection.method = 'vst', nfeatures = 2000)
#top10_varfeat <- head(VariableFeatures(st_hcc05), 10)
# "IGKC"       "MALAT1"     "IGLC2"      "SLA"        "LINC01128"  "SLC5A4-AS1" "IGHG1"      "TBX15" "CHML"       "SPART" 

st_hcc05 <- FindSpatiallyVariableFeatures(st_hcc05, assay = "SCT", features = VariableFeatures(st_hcc05)[1:500], 
                                          selection.method = "markvariogram")
#top_features <- head(SpatiallyVariableFeatures(st_hcc05, selection.method = "markvariogram"), 10)
# "MALAT1"     "CYP2E1"     "HP"         "MT-CO1"     "ALKBH3"     "CD96"       "STC2"       "AL355490.1" "TANC2"      "C20orf194"

##### SCALING THE DATA
all_genes <- rownames(st_hcc05)
st_hcc05 <- ScaleData(st_hcc05, features = all_genes)

##### DIMENSIONALITY REDUCTION, CLUSTERING, VISUALIZATION
st_hcc05 <- RunPCA(st_hcc05, assay = "SCT", verbose = F)
print(st_hcc05[['pca']])

# run UMAP for dimensional reduction
st_hcc05 <- RunUMAP(st_hcc05, reduction = "pca", dims = 1:50)
DimPlot(st_hcc05, reduction = "umap", label = T, split.by = 'new.ident')

# This whole tissue slide is tumor region according to the pathological analysis
SpatialDimPlot(st_hcc05, cells.highlight = CellsByIdentities(object = st_hcc05, idents = c(0, 1, 2, 3)), 
               facet.highlight = T, ncol = 2)

# Find highly expressed genes in each cluster
hcc05_markers <- FindAllMarkers(st_hcc05, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)

markers_df <- as.data.frame(hcc05_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 10, order_by = avg_log2FC))

de_markers <- FindMarkers(st_hcc05, ident.1 = 0, ident.2 = 2)
de_markers <- de_markers[order(de_markers$avg_log2FC),]

SpatialFeaturePlot(object = st_hcc05, features = c("MALAT1",'NDUFB4','CLDN2'), 
                   alpha = c(0.1, 1), ncol = 3)

write.csv(markers_df,
          sprintf("%s/hcc05_markers_stlearn.csv",output_dir),
          row.names = FALSE)


current.cluster.ids <- c(0, 1, 2, 3)
new.cluster.ids <- c("Tumor-1", "Tumor-2 (MALAT1)", "Tumor-3 (NDUFB4)", "Tumor-4 (CLDN2)")
st_hcc05@meta.data$seurat_clusters <- plyr::mapvalues(x = st_hcc05@meta.data$new.ident, from = current.cluster.ids, to = new.cluster.ids)
annotate_meta <- st_hcc05@meta.data

st_hcc05@meta.data <- annotate_meta

SpatialDimPlot(st_hcc05, group.by = "seurat_clusters") + 
  theme(text = element_text(size = 20)) &
  guides(fill=guide_legend(override.aes=list(size=5)))

DoHeatmap(object = st_hcc05, 
          features = markers_df$gene,
          group.by = 'seurat_clusters')  + 
  theme(text = element_text(size = 20))


write.csv(st_hcc05@meta.data,
          sprintf("%s/st_hcc05_metadata.csv",output_dir),
          row.names = TRUE)

save(st_hcc05, 
     file = sprintf("%s/hcc05_UMAPreduced.Rda",output_dir))
######################  Result plotting ####################
SpatialDimPlot(st_hcc05, group.by = 'seurat_clusters')

auc_mtx <- (read.table(sprintf('./%s/auc_mtx.csv', SCENIC_output_dir), header = TRUE, row.names = 1, 
                      stringsAsFactors = FALSE, sep = ',')) #load pyscenic output

#subsetting important transcription factor modules for downstream analysis
auc_sub <- auc_mtx[ c('ETS1_', 'STAT3_', 'MYC_', 'HEY2_', 'ETS2_', 'FLI1_', 'HEY1_', 'EGR1_', 'NR2F2_', 'JUN_', 'FOS_', 'XBP1_'), ]
auc_sub <- t(auc_sub)
colnames(auc_sub) <- c('ETS1_Activity', 'STAT3_Activity', 'MYC_Activity', 'HEY2_Activity',
                       'ETS2_Activity', 'FLI1_Activity', 'HEY1_Activity', 'EGR1_Activity',
                       'NR2F2_Activity', 'JUN_Activity', 'FOS_Activity', 'XBP1_Activity')

#Merge transcription factor module activities into metadata of the seurat for plotting
st_hcc05@meta.data <- merge(st_hcc05@meta.data, auc_sub,
                            by = 'row.names', all = TRUE)

rownames(st_hcc05@meta.data) <- st_hcc05@meta.data$Row.names

#Remove spot with extreme values
st_hcc05 <- subset(x = st_hcc05, subset = MYC_Activity < 0.05)
st_hcc05 <- subset(x = st_hcc05, subset = STAT3_Activity < 0.3)
st_hcc05 <- subset(x = st_hcc05, subset = ETS1_Activity < 0.2)
st_hcc05 <- subset(x = st_hcc05, subset = HEY2_Activity < 0.15)
st_hcc05 <- subset(x = st_hcc05, subset = FLI1_Activity < 0.2)

#Spatial plot of transcription factor activities.
SpatialFeaturePlot(object = st_hcc05, 
                   features = c('ETS1_Activity','ETS2_Activity', 'STAT3_Activity',
                                'MYC_Activity', 'HEY2_Activity', 'FLI1_Activity') , 
                   alpha = c(0.1, 1), ncol=2) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")

SpatialFeaturePlot(object = st_hcc05, 
                   features = c('ETS1_Activity','ETS2_Activity', 'STAT3_Activity',
                                'MYC_Activity', 'HEY2_Activity', 'FLI1_Activity') , 
                   alpha = c(0.1, 1), ncol=3) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


SpatialFeaturePlot(object = st_hcc05, features = c('JAG1', 'DLL1', 'NOTCH3', 'NOTCH4'), 
                   alpha = c(0.1, 1), ncol = 2) & scale_fill_continuous(low = "white", high = "red") & theme(legend.position='right') 



#-- Module Scores ---------------------------------------------
CAF_list <- list(c('COL6A3', 'COL3A1', 'COL1A1', 'COL1A2' ,'VIM', 'DCN'))
st_hcc05 <- AddModuleScore(object = st_hcc05, features = CAF_list, name = "CAF_score")
SpatialFeaturePlot(object = st_hcc05, features = "CAF_score1", alpha = c(0.1, 1))

BCellList <- list(c('CD19', 'CD79A', 'CD40', 'CR2', 'FCER2A',
                    'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                    'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'))
st_hcc05 <- AddModuleScore(object = st_hcc05, features = BCellList, name = "Bcell_score")
SpatialFeaturePlot(object = st_hcc05, features = "Bcell_score1", alpha = c(0.1, 1))


plasemaCellList <- list(c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                          'IGHM', 'IGHD', 'IGHG3','IGHG1',
                          'IGKC','IGLC1','IGLC2', 'IGHA1'))
st_hcc05 <- AddModuleScore(object = st_hcc05, features = plasemaCellList, name = "Plasma_score")
SpatialFeaturePlot(object = st_hcc05, features = "Plasma_score1", alpha = c(0.1, 1))

angioList <- list(c('PECAM1', 'VWF', 'ENG', 'CDH5', 'CD34'))
st_hcc05 <- AddModuleScore(object = st_hcc05, features = angioList, name = "VasScore")
SpatialFeaturePlot(object = st_hcc05, features = "VasScore1", alpha = c(0.1, 1))

#-- Spatial Plots ---------------------------------------------
SpatialFeaturePlot(object = st_hcc05, features = c('ANG',  'FGA', 'HP', 'C3',
                                                   'APOA1', 'AFP', 'VEGFA',
                                                   'ALDH2', 'TGFB1', 'IL6R'), 
                   alpha = c(0.1, 1), ncol = 4) & scale_fill_continuous(low = "white", high = "red")


SpatialFeaturePlot(object = st_hcc05, features = c('ALDH1A1', 'ALDH1A3', 'ALDH1A2', 'ALDH2'), 
                   alpha = c(0.1, 1), ncol = 4) & scale_fill_continuous(low = "white", high = "red")

SpatialFeaturePlot(object = st_hcc05, 
                   features = c('COL1A1', 'COL1A2', 'COL6A3', 'COL6A1', 'COL3A1','VIM', 'THY1', 'DCN') , 
                   alpha = c(0.1, 1), ncol=3) & scale_fill_continuous(low = "white", high = "red")



########################## Angiokine ###########################

angiokine <- c('ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
               'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
               'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
               'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')

angiokine_hcc05 <- c( 'ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
                     'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
                     'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
                     'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')


VlnPlot(object = st_hcc05, features = angiokine,
        split.by = 'seurat_clusters', ncol=6, pt.size=0) & 
  stat_compare_means(size=4)  & 
  theme(text = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')

