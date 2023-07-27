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
raw_data_path <- "./PATH_TO_HCC-7_RAW"
ST_result_path <- "./PATH_TO_STLEARN_RESULT"
SCENIC_output_dir <- './PATH_TO_SCENIC_OUTPUT/'
setwd(work_dir)

st_hcc07 <- Load10X_Spatial(raw_data_path)
data_SME_HCC07_identity <- read.csv(sprintf("%s/data_SME_hcc07_identity.csv", ST_result_path), header = T)

######################## Pre-processing #########################
new_meta <- st_hcc07@meta.data
new_meta$new.ident <- data_SME_HCC07_identity$X_pca_kmeans[which(rownames(new_meta) %in% data_SME_HCC07_identity$X)]
st_hcc07@meta.data <- new_meta

## Visualize clustering
Idents(st_hcc07) <- 'new.ident'
SpatialDimPlot(st_hcc07, cells.highlight = CellsByIdentities(object = st_hcc07, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), 
               facet.highlight = T, ncol = 4)

## QC and filtering
st_hcc07[['percent_mt']] <- PercentageFeatureSet(st_hcc07, pattern = '^MT-')
st_hcc07 <- subset(st_hcc07, subset = nFeature_Spatial > 100 & percent_mt < 40)
# removing spots with zero counts
st_hcc07 <- subset(st_hcc07, subset = nCount_Spatial > 0)
# normalizing with SCTransform
st_hcc07 <- SCTransform(st_hcc07, assay = "Spatial", vars.to.regress = 'percent_mt',verbose = F)

# FindSpatiallyVariables search for features without pre-annotation
# "markvariogram" finds genes whose expression depends on spatial location
st_hcc07 <- FindVariableFeatures(st_hcc07, selection.method = 'vst', nfeatures = 2000)
#top10_varfeat <- head(VariableFeatures(st_hcc07), 10)
# "IGKC"   "IGHG3"  "COL1A1" "COL1A2" "HBB"    "IGHA1"  "COL3A1" "IGFBP1" "IGLC1"  "HBA2"  

st_hcc07 <- FindSpatiallyVariableFeatures(st_hcc07, assay = "SCT", features = VariableFeatures(st_hcc07)[1:500], 
                                          selection.method = "markvariogram")
#top_features <- head(SpatiallyVariableFeatures(st_hcc07, selection.method = "markvariogram"), 10)
# "MT-CO1"   "COL1A2"   "COL1A1"   "APOA2"    "MALAT1"   "SERPINA1" "COL3A1"   "SPINK1"   "AFP"      "MGP"     

##### SCALING THE DATA
all_genes <- rownames(st_hcc07)
st_hcc07 <- ScaleData(st_hcc07, features = all_genes)

##### DIMENSIONALITY REDUCTION, CLUSTERING, VISUALIZATION
st_hcc07 <- RunPCA(st_hcc07, assay = "SCT", verbose = F)
print(st_hcc07[['pca']])

# run UMAP for dimensional reduction
st_hcc07 <- RunUMAP(st_hcc07, reduction = "pca", dims = 1:50)
DimPlot(st_hcc07, reduction = "umap", label = T, split.by = 'new.ident')

# Subset out cluster0 4 6 7 (little spots)
st_hcc07 <- subset(st_hcc07, idents = c(1, 2, 3, 5))
SpatialDimPlot(st_hcc07, cells.highlight = CellsByIdentities(object = st_hcc07, idents = c(1, 2, 3, 5)), 
               facet.highlight = T, ncol = 2)

# Find highly expressed genes in each cluster
hcc07_markers <- FindAllMarkers(st_hcc07, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
markers_df <- as.data.frame(hcc07_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 10, order_by = avg_log2FC))
write.csv(markers_df,
          sprintf("%s/hcc07_markers_stlearn.csv",output_dir),
          row.names = FALSE)

#cluster 5 AFP AC132217.1 RELN ENHO RPS23 - 
#cluster 2 PLA2G2A SLPI PGC SPINK1 INHA - CAFs
#cluster 1 COL1A1 COL1A2 COL3A1 SPARC TIMP1 - 
#cluster 3 MALAT1 MT-ND2 MT-ND3 MT-ATP6 NEAT1 - 

### Annotate each cell cluster
current.cluster.ids <- c(1,2,3,5)
new.cluster.ids <- c("Tumor-1 (PLA2G2A)", "CAFs", "Tumor-2 (MALAT1)", "Tumor-3 (AFP)")
st_hcc07@meta.data$seurat_clusters <- plyr::mapvalues(x = st_hcc07@meta.data$new.ident, from = current.cluster.ids, to = new.cluster.ids)

SpatialDimPlot(st_hcc07, group.by = "seurat_clusters") + 
  theme(text = element_text(size = 20)) &
  guides(fill=guide_legend(override.aes=list(size=5)))

DoHeatmap(object = st_hcc07, 
          features = markers_df$gene,
          group.by = 'seurat_clusters')  + 
  theme(text = element_text(size = 20))

write.csv(st_hcc07@meta.data,
          sprintf("%s/st_hcc07_metadata.csv",output_dir),
          row.names = TRUE)

SpatialFeaturePlot(object = st_hcc07, features = c('AFP','COL1A1','PLA2G2A','MALAT1','ASS1'), 
                   alpha = c(0.1, 1), ncol = 3)

annotate_meta <- st_hcc07@meta.data
annotate_meta$annotate.ident <- combined_meta$seurat_clusters[which(rownames(annotate_meta) %in% rownames(combined_meta))]
st_hcc07@meta.data <- annotate_meta
SpatialDimPlot(st_hcc07, group.by = 'annotate.ident')
SpatialDimPlot(st_hcc07, group.by = 'seurat_clusters')

save(st_hcc07, 
     file = sprintf("%s/hcc07_UMAPreduced.Rda",output_dir))


############################   TF module result Plotting  ################################

SpatialDimPlot(st_hcc07, group.by = "seurat_clusters")

auc_mtx = (read.table(sprintf('./%s/auc_mtx.csv', SCENIC_output_dir), header = TRUE, row.names = 1, 
                      stringsAsFactors = FALSE, sep = ',')) #load pyscenic output


auc_sub <- auc_mtx[ c('ETS1_', 'STAT3_', 'MYC_',
                      'HEY2_', 'ETS2_', 'FLI1_',
                      'EGR1_', 'NR2F2_'), ]

auc_sub <- t(auc_sub)
colnames(auc_sub) <- c('ETS1_Activity', 'STAT3_Activity', 'MYC_Activity', 
                       'HEY2_Activity', 'ETS2_Activity', 'FLI1_Activity',
                       'EGR1_Activity', 'NR2F2_Activity')

st_hcc07@meta.data <- merge(st_hcc07@meta.data, auc_sub,
                            by = 'row.names', all = TRUE)

rownames(st_hcc07@meta.data) <- st_hcc07@meta.data$Row.names

#------------------------      subsetting by cebpb score ----------------------
st_hcc07 <- subset(x = st_hcc07, subset = MYC_Activity < 0.05)
st_hcc07 <- subset(x = st_hcc07, subset = HEY2_Activity < 0.3)
st_hcc07 <- subset(x = st_hcc07, subset = STAT3_Activity < 0.2)
st_hcc07 <- subset(x = st_hcc07, subset = ETS2_Activity < 0.2)
st_hcc07 <- subset(x = st_hcc07, subset = FLI1_Activity < 0.08)

##############################################################################################################
##HCC09 Figure D Spatially resolved transcription factor activities of ETS1, MYC, HEY2, FLI1, ETS2, and ETS1.#
##############################################################################################################
SpatialFeaturePlot(object = st_hcc07, 
                   features = c('ETS1_Activity','ETS2_Activity', 'STAT3_Activity',
                                'MYC_Activity', 'HEY2_Activity', 'FLI1_Activity') , 
                   alpha = c(0.1, 1), ncol=3) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


VlnPlot(object = st_hcc07, features = c('CEBPB_Activity','JUN_Activity', 'MYC_Activity', 'ETS1_Activity'),
        split.by = 'seurat_clusters', ncol=2, pt.size=0) & 
  stat_compare_means(size=5)  & 
  theme(text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')

#-- Module Scores ---------------------------------------------
CAF_list <- list(c('COL6A3', 'COL3A1', 'COL1A1', 'COL1A2' ,'VIM', 'DCN'))
st_hcc07 <- AddModuleScore(object = st_hcc07, features = CAF_list, name = "CAF_score")
SpatialFeaturePlot(object = st_hcc07, features = "CAF_score1", alpha = c(0.1, 1))

BCellList <- list(c('CD19', 'CD79A', 'CD40', 'CR2', 'FCER2A',
                    'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                    'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'))
st_hcc07 <- AddModuleScore(object = st_hcc07, features = BCellList, name = "Bcell_score")
SpatialFeaturePlot(object = st_hcc07, features = "Bcell_score1", alpha = c(0.1, 1))


plasemaCellList <- list(c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                          'IGHM', 'IGHD', 'IGHG3','IGHG1',
                          'IGKC','IGLC1','IGLC2', 'IGHA1'))
st_hcc07 <- AddModuleScore(object = st_hcc07, features = plasemaCellList, name = "Plasma_score")
SpatialFeaturePlot(object = st_hcc07, features = "Plasma_score1", alpha = c(0.1, 1))


angioList <- list(c('PECAM1', 'VWF', 'ENG', 'CDH5', 'CD34'))
st_hcc07 <- AddModuleScore(object = st_hcc07, features = angioList, name = "VasScore")
SpatialFeaturePlot(object = st_hcc07, features = "VasScore1", alpha = c(0.1, 1))

#-- Spatial Plots ---------------------------------------------
SpatialFeaturePlot(object = st_hcc07, features = c('ANG',  'FGA', 'HP', 'C3',
                                                   'APOA1', 'AFP', 'VEGFA',
                                                   'ALDH2', 'TGFB1', 'IL6R'), 
                   alpha = c(0.1, 1), ncol = 4) & scale_fill_continuous(low = "white", high = "red")


SpatialFeaturePlot(object = st_hcc07, features = c('ALDH1A1', 'ALDH1A3', 'ALDH1A2', 'ALDH2'), 
                   alpha = c(0.1, 1), ncol = 4) & scale_fill_continuous(low = "white", high = "red")

SpatialFeaturePlot(object = st_hcc07, 
                   features = c('COL1A1', 'COL1A2', 'COL6A3', 'COL6A1', 'COL3A1','VIM', 'THY1', 'DCN') , 
                   alpha = c(0.1, 1), ncol=3) & scale_fill_continuous(low = "white", high = "red")


SpatialFeaturePlot(object = st_hcc07, 
                   features = c('ETS1_Activity','ETS2_Activity', 'STAT3_Activity',
                                'MYC_Activity', 'HEY2_Activity', 'FLI1_Activity') , 
                   alpha = c(0.1, 1), ncol=3) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


################### subsetting ###########
st_hcc07_tumor <- subset(x = st_hcc07,  idents = c(1, 3, 5))
st_hcc07_sub <- subset(x = st_hcc07,  idents = c(2, 3))


scatterPlot(st_hcc07_sub, 'NR2F2', 'VasScore1', 'Vas Score') 


scatterPlot(st_hcc07_sub, 'EGR1', 'VasScore1', 'Vas Score') 

SpatialFeaturePlot(object = st_hcc07, features = c('JAG1', 'DLL1', 'NOTCH3', 'NOTCH4'), 
                   alpha = c(0.1, 1), ncol = 2) & scale_fill_continuous(low = "white", high = "red") & theme(legend.position='right') 



########################## Angiokine ###########################
angiokine <- c('ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
               'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
               'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
               'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')

angiokine_hcc07 <- c( 'MET',  'NRAS', 'VEGFA', 'VEGFB')


VlnPlot(object = st_hcc07, features = angiokine,
        split.by = 'seurat_clusters', ncol=6, pt.size=0) & 
  stat_compare_means(size=4)  & 
  theme(text = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')


#####################################          GESA              ###########################

# Find highly expressed genes in each cluster
hcc09_markers <- FindAllMarkers(st_hcc07, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
markers_df <- as.data.frame(hcc09_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 5, order_by = avg_log2FC))

library(clusterProfiler)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(AnnotationDbi)
# Extract the foldchanges for GSEA
res <- hcc09_markers[hcc09_markers$cluster == 2,]
foldchanges <- res$avg_log2FC
names(foldchanges) <- res$gene
#names(foldchanges) <- rownames(res)
foldchanges <- sort(foldchanges, decreasing = TRUE)
# use mapIds method to obtain Entrez IDs
names(foldchanges) <- mapIds(org.Hs.eg.db, names(foldchanges), 'ENTREZID', 'SYMBOL')
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 1000, # default number permutations
                    minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

