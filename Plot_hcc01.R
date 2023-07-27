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
raw_data_path <- "./PATH_TO_HCC-1_RAW"
ST_result_path <- "./PATH_TO_STLEARN_RESULT"
SCENIC_output_dir <- './PATH_TO_SCENIC_OUTPUT/'
setwd(work_dir)

st_hcc01 <- Load10X_Spatial(raw_data_path)
data_SME_HCC01_identity <- read.csv(sprintf("%s/data_SME_hcc01_identity.csv", ST_result_path), header = T)

######################## Pre-processing #########################
# Assign cell population info generated with StLearn
new_meta <- st_hcc01@meta.data
new_meta$new.ident <- data_SME_HCC01_identity$X_pca_kmeans[which(rownames(new_meta) %in% data_SME_HCC01_identity$X)]
st_hcc01@meta.data <- new_meta

## Visualize clustering
Idents(st_hcc01) <- 'new.ident'
SpatialDimPlot(st_hcc01, cells.highlight = CellsByIdentities(object = st_hcc01, idents = c(0, 1, 2, 3, 4, 5,6,7,8)), 
               facet.highlight = T, ncol = 3)

## QC and filtering
st_hcc01[['percent_mt']] <- PercentageFeatureSet(st_hcc01, pattern = '^MT-')
st_hcc01 <- subset(st_hcc01, subset = nFeature_Spatial > 100 & percent_mt < 40)
# removing spots with zero counts
st_hcc01 <- subset(st_hcc01, subset = nCount_Spatial > 0)
# normalizing with SCTransform
st_hcc01 <- SCTransform(st_hcc01, assay = "Spatial", vars.to.regress = 'percent_mt',verbose = F)

# FindSpatiallyVariables search for features without pre-annotation
# "markvariogram" finds genes whose expression depends on spatial location
st_hcc01 <- FindVariableFeatures(st_hcc01, selection.method = 'vst', nfeatures = 2000)

st_hcc01 <- FindSpatiallyVariableFeatures(st_hcc01, assay = "SCT", features = VariableFeatures(st_hcc01)[1:500], 
                                          selection.method = "markvariogram")

##### SCALING THE DATA
all_genes <- rownames(st_hcc01)
st_hcc01 <- ScaleData(st_hcc01, features = all_genes)

##### DIMENSIONALITY REDUCTION, CLUSTERING, VISUALIZATION
st_hcc01 <- RunPCA(st_hcc01, assay = "SCT", verbose = F)
print(st_hcc01[['pca']])

# run UMAP for dimensional reduction
st_hcc01 <- RunUMAP(st_hcc01, reduction = "pca", dims = 1:50)
DimPlot(st_hcc01, reduction = "umap", label = T, split.by = 'new.ident')

# Subset out cluster678 (little spots)
st_hcc01 <- subset(st_hcc01, idents = c(0, 1, 2, 3, 4,5))
SpatialDimPlot(st_hcc01, cells.highlight = CellsByIdentities(object = st_hcc01, idents = c(0, 1, 2, 3, 4,5)), 
               facet.highlight = T, ncol = 3)

# Find highly expressed genes in each cluster
hcc01_markers <- FindAllMarkers(st_hcc01, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
markers_df <- as.data.frame(hcc01_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 10, order_by = avg_log2FC))
write.csv(markers_df,
          sprintf("%s/hcc01_markers_stlearn.csv",output_dir),
          row.names = FALSE)

#1 IGLC1 IGLC2 IGKC IGLC3 IGHG3 - B
#5 GGH AKR1C1 SPP1 AKR1B10 UGT2B4 - immune+tumor
#2 AKR1C4 ITIH2 NDUFB11 MAT1A C2 - DCs (NDUFB11,TYROBP,HLA-DPA1,CUTA)
#4 APOD MGP SFRP2 COL1A1 IGFBP4 - CAFs
#3 IGLC2 IGKC IGLC3 CD74 TMSB4X - immune
#0 ALB CRP APOA1 TF SAA1 - tumor

### Annotate each cell cluster
current.cluster.ids <- c(0, 1, 2, 3, 4,5)
new.cluster.ids <- c("Tumor-1 (ALB)", "B cell", "DCs", "Immune", "CAFs", 'Tumor-2 (GGH)')
st_hcc01@meta.data$seurat_clusters <- plyr::mapvalues(x = st_hcc01@meta.data$new.ident, from = current.cluster.ids, to = new.cluster.ids)

SpatialDimPlot(st_hcc01, group.by = "seurat_clusters") + 
  theme(text = element_text(size = 20)) &
  guides(fill=guide_legend(override.aes=list(size=5)))

DoHeatmap(object = st_hcc01, 
          features = markers_df$gene,
          group.by = 'seurat_clusters')  + 
  theme(text = element_text(size = 20))


write.csv(st_hcc01@meta.data,
          sprintf("%s/st_hcc01_metadata.csv",output_dir),
          row.names = TRUE)

VlnPlot(st_hcc01, features = c("ITGAX",'NDUFB11','TYROBP','HLA-DPA1','CUTA'), split.by = 'new.ident')
VlnPlot(st_hcc01, features = c("IGLC1",'IGLC2','IGKC','IGLC3','IGHG3'), split.by = 'new.ident')

SpatialFeaturePlot(object = st_hcc01, features = c("ITGAX",'NDUFB11','TYROBP','HLA-DPA1','CUTA'), 
                   alpha = c(0.1, 1), ncol = 3)


save(st_hcc01, 
     file = sprintf("%s/hcc01_UMAPreduced.Rda",output_dir))
## HCC01plot
#########################           plot           ###############################

SpatialDimPlot(st_hcc01, group.by = "seurat_clusters")

auc_mtx = (read.table(sprintf('./%s/auc_mtx.csv', SCENIC_output_dir), header = TRUE, row.names = 1, 
                      stringsAsFactors = FALSE, sep = ',')) #load pyscenic output

auc_sub <- auc_mtx[c('PAX5_', 'TCF4_', 'ETS1_', 
                     'IRF4_', 'XBP1_', 'IRF1_',
                     'JUN_', 'FOS_', 'EGR1_', 
                     'POU5F1_', 'NR2F2_'), ]
auc_sub <- t(auc_sub)
colnames(auc_sub) <- c('PAX5_Activity', 'TCF4_Activity', 'ETS1_Activity',
                       'IRF4_Activity', 'XBP1_Activity', 'IRF1_Activity',
                       'JUN_Activity', 'FOS_Activity', 'EGR1_Activity', 'POU5F1_Activity', 'NR2F2_Activity')


st_hcc01@meta.data <- merge(st_hcc01@meta.data, auc_sub,
                            by = 'row.names', all = TRUE)

rownames(st_hcc01@meta.data) <- st_hcc01@meta.data$Row.names
#-- Module Scores ---------------------------------------------
CAF_list <- list(c('COL6A3', 'COL3A1', 'COL1A1', 'COL1A2' ,'VIM', 'DCN'))
st_hcc01 <- AddModuleScore(object = st_hcc01, features = CAF_list, name = "CAF_score")
SpatialFeaturePlot(object = st_hcc01, features = "CAF_score1", alpha = c(0.1, 1))

#CD21: CR2
BCellList <- list(c('CD19', 'CD79A', 'CD40', 'CR2',
                    'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                    'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'))
st_hcc01 <- AddModuleScore(object = st_hcc01, features = BCellList, name = "Bcell_score")
SpatialFeaturePlot(object = st_hcc01, features = "Bcell_score1", alpha = c(0.1, 1))

plasemaCellList <- list(c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                          'IGHM', 'IGHD', 'IGHG3','IGHG1',
                          'IGKC','IGLC1','IGLC2', 'IGHA1'))
st_hcc01 <- AddModuleScore(object = st_hcc01, features = plasemaCellList, name = "Plasma_score")
SpatialFeaturePlot(object = st_hcc01, features = "Plasma_score1", alpha = c(0.1, 1))

angioList <- list(c('PECAM1', 'VWF', 'ENG', 'CDH5', 'CD34'))
st_hcc01 <- AddModuleScore(object = st_hcc01, features = angioList, name = "VasScore")
SpatialFeaturePlot(object = st_hcc01, features = "VasScore1", alpha = c(0.1, 1))

SpatialDimPlot(st_hcc01, group.by = 'seurat_clusters')

############################# subsetting ####################################
#Subsetting out extreme values
st_hcc01 <- subset(x = st_hcc01, subset = IRF1_Activity < 0.07)
st_hcc01 <- subset(x = st_hcc01, subset = JUN_Activity < 0.07)
st_hcc01 <- subset(x = st_hcc01, subset = PAX5_Activity < 0.1)

################################          Spatial           ########################################


SpatialFeaturePlot(object = st_hcc01, 
                   features = c('CD79A', 'CD19', 'COL6A3', 'COL3A1', 'COL1A1',
                                'CD22', 'CIITA','COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=5) & 
  theme(legend.position='right')


SpatialFeaturePlot(object = tumor_region, 
                   features = c('IGHG3', 'IGHG1', 'IGLC1', 'IGLC2',
                                'CCL21', 'CCL19') , 
                   alpha = c(0.1, 1), ncol=4) & scale_fill_continuous(low = "white", high = "blue")






#-- spatial TF HCC01 Figure A-------------

SpatialFeaturePlot(object = st_hcc01, 
                   features = c('JUN_Activity') , 
                   alpha = c(0.1, 1)) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


SpatialFeaturePlot(object = st_hcc01, 
                   features = c('FOS_Activity') , 
                   alpha = c(0.1, 1)) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


SpatialFeaturePlot(object = st_hcc01, 
                   features = c('PAX5_Activity') , 
                   alpha = c(0.1, 1)) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")




SpatialFeaturePlot(object = st_hcc01, 
                   features = c('CD79A', 'CD19','CD22', 'CIITA') , 
                   alpha = c(0.2, 1), ncol=2) & 
  theme(legend.position='right') 


SpatialFeaturePlot(object = st_hcc01, 
                   features = c('COL6A3', 'COL3A1', 'COL1A1',
                                'COL1A2' ,'VIM', 'DCN') , 
                   alpha = c(0.2, 1), ncol=3) & 
  theme(legend.position='right') 



tumor_region <-  subset(st_hcc01, idents = c(0, 5))


SpatialFeaturePlot(object = tumor_region, 
                   features = c('VEGFA', 'VEGFB', 'ANG',
                                'ANGPT2', 'CD34',  'PECAM1') , 
                   alpha = c(0.2, 1), ncol=3) & 
  theme(legend.position='right')

SpatialFeaturePlot(object = st_hcc01, 
                   features = c('HAMP', 'CEBPA', 'HIF1A',
                                'CD24') , 
                   alpha = c(0.2, 1), ncol=3) & 
  theme(legend.position='right')


VlnPlot(object = tumor_region, features = c('HAMP', 'CEBPA', 'HIF1A',
                                            'CD24'),
        split.by = 'seurat_clusters', ncol=2) & 
  stat_compare_means(size=5)  & 
  theme(text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')

###########################################################################################
#Supplemental Figure 2 B
st_hcc01_stroma <-  subset(st_hcc01, idents = c(3, 4))

violin_list = FetchData(object = st_hcc01_stroma, vars = c('VasScore1', 'CAF_score1', 'Bcell_score1', 'Plasma_score1', 'seurat_clusters'))

my_comparisons <- list( c("CAFs", "Immune") )

p_vas <- ggplot(violin_list, aes(x= seurat_clusters, y = VasScore1)) + 
  geom_violin(aes(fill=seurat_clusters)) +     
  theme_minimal() +
  geom_boxplot(width=0.1) + labs(y= "Vaculature Score", x = "", fill='Cluster') +
  guides(col=guide_legend("Cluster")) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 16, hjust= 0.5, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 16, color = 'black'),
        legend.text = element_text(size = 14)) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 5)

p_caf <- ggplot(violin_list, aes(x= seurat_clusters, y = CAF_score1)) + 
  geom_violin(aes(fill=seurat_clusters)) +     
  theme_minimal() +
  geom_boxplot(width=0.1) + labs(y= "ECM Score", x = "", fill='Cluster') +
  guides(col=guide_legend("Cluster")) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 16, hjust= 0.5, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 16, color = 'black'),
        legend.text = element_text(size = 14)) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 5)

p_bcell <- ggplot(violin_list, aes(x= seurat_clusters, y = Bcell_score1)) + 
  geom_violin(aes(fill=seurat_clusters)) +     
  theme_minimal() +
  geom_boxplot(width=0.1) + labs(y= "B cell Score", x = "", fill='Cluster') +
  guides(col=guide_legend("Cluster")) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 16, hjust= 0.5, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 16, color = 'black'),
        legend.text = element_text(size = 14)) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 5)

p_plasma <- ggplot(violin_list, aes(x= seurat_clusters, y = Plasma_score1)) + 
  geom_violin(aes(fill=seurat_clusters)) +     
  theme_minimal() +
  geom_boxplot(width=0.1) + labs(y= "Plasma Score", x = "", fill='Cluster') +
  guides(col=guide_legend("Cluster")) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 16, hjust= 0.5, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 16, color = 'black'),
        legend.text = element_text(size = 14)) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 5)



p_bcell + p_plasma + p_caf + p_vas






#####################################          GESA              ###########################

library(clusterProfiler)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(AnnotationDbi)
# Extract the foldchanges for GSEA
res <- hcc01_markers[hcc01_markers$cluster == 3,]
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
