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
raw_data_path <- "./PATH_TO_HCC-6_RAW"
ST_result_path <- "./PATH_TO_STLEARN_RESULT"
SCENIC_output_dir <- './PATH_TO_SCENIC_OUTPUT/'
setwd(work_dir)

st_hcc06 <- Load10X_Spatial(raw_data_path)
data_SME_HCC06_identity <- read.csv(sprintf("%s/data_SME_hcc06_identity.csv", ST_result_path), header = T)

######################## Pre-processing #########################
# Assign cell population info generated with StLearn
new_meta <- st_hcc06@meta.data
new_meta$new.ident <- data_SME_HCC06_identity$X_pca_kmeans[which(rownames(new_meta) %in% data_SME_HCC06_identity$X)]
st_hcc06@meta.data <- new_meta

## Visualize clustering
Idents(st_hcc06) <- 'new.ident'
SpatialDimPlot(st_hcc06, cells.highlight = CellsByIdentities(object = st_hcc06, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8)), 
               facet.highlight = T, ncol = 3)

## QC and filtering
st_hcc06[['percent_mt']] <- PercentageFeatureSet(st_hcc06, pattern = '^MT-')
st_hcc06 <- subset(st_hcc06, subset = nFeature_Spatial > 100 & percent_mt < 40)
# removing spots with zero counts
st_hcc06 <- subset(st_hcc06, subset = nCount_Spatial > 0)
# normalizing with SCTransform
st_hcc06 <- SCTransform(st_hcc06, assay = "Spatial", vars.to.regress = 'percent_mt',verbose = F)

# FindSpatiallyVariables search for features without pre-annotation
# "markvariogram" finds genes whose expression depends on spatial location
st_hcc06 <- FindVariableFeatures(st_hcc06, selection.method = 'vst', nfeatures = 2000)
#top10_varfeat <- head(VariableFeatures(st_hcc06), 10)
# "IGKC"   "ALB"    "IGHG1"  "TIMP1"  "COL1A1" "IGLC2"  "IGHG3"  "H19"    "APOC3"  "HBA2"  

st_hcc06 <- FindSpatiallyVariableFeatures(st_hcc06, assay = "SCT", features = VariableFeatures(st_hcc06)[1:500], 
                                          selection.method = "markvariogram")
#top_features <- head(SpatiallyVariableFeatures(st_hcc06, selection.method = "markvariogram"), 10)
# "MT-ND2"   "MT-CO1"   "TF"       "TIMP1"    "ALB"      "DLK1"     "COL1A1"   "FTL"      "SERPINA1" "APOC3" 

##### SCALING THE DATA
all_genes <- rownames(st_hcc06)
st_hcc06 <- ScaleData(st_hcc06, features = all_genes)

##### DIMENSIONALITY REDUCTION, CLUSTERING, VISUALIZATION
st_hcc06 <- RunPCA(st_hcc06, assay = "SCT", verbose = F)
print(st_hcc06[['pca']])

# run UMAP for dimensional reduction
st_hcc06 <- RunUMAP(st_hcc06, reduction = "pca", dims = 1:50)
DimPlot(st_hcc06, reduction = "umap", label = T)

# Subset out cluster5 6 7 (little spots)
st_hcc06 <- subset(st_hcc06, idents = c(0, 1, 2, 4, 8))
SpatialDimPlot(st_hcc06, cells.highlight = CellsByIdentities(object = st_hcc06, idents = c(0, 1, 2, 4, 8)), 
               facet.highlight = T, ncol = 3)

# Find highly expressed genes in each cluster
hcc06_markers <- FindAllMarkers(st_hcc06, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
markers_df <- as.data.frame(hcc06_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 10, order_by = avg_log2FC))
write.csv(markers_df,
          sprintf("%s/hcc06_markers_stlearn.csv",output_dir),
          row.names = FALSE)

### Annotate each cell cluster
current.cluster.ids <- c(0, 1, 2, 4, 8)
new.cluster.ids <- c("Tumor-1 (LRRC75A)", "Tumor-2 (PLA2G1B)", "Tumor-3 (CYP2E1)", "CAFs", 'Immune')
st_hcc06@meta.data$seurat_clusters <- plyr::mapvalues(x = st_hcc06@meta.data$new.ident, from = current.cluster.ids, to = new.cluster.ids)


SpatialFeaturePlot(object = st_hcc06, features = c('ALB','LRRC75A','PLA2G1B','IGKC','COL1A1'), 
                   alpha = c(0.1, 1), ncol = 3)

SpatialDimPlot(st_hcc06, group.by = "seurat_clusters") + 
  theme(text = element_text(size = 20)) &
  guides(fill=guide_legend(override.aes=list(size=5)))

write.csv(st_hcc06@meta.data,
          sprintf("%s/st_hcc06_metadata.csv",output_dir),
          row.names = TRUE)

DoHeatmap(object = st_hcc06, 
          features = markers_df$gene,
          group.by = 'seurat_clusters')  + 
  theme(text = element_text(size = 20))

save(st_hcc06, 
     file = sprintf("%s/HCC06_UMAPreduced.Rda",output_dir))
############################### Manuscript #########################################
auc_mtx = (read.table(sprintf('./%s/auc_mtx.csv', SCENIC_output_dir), header = TRUE, row.names = 1, 
                      stringsAsFactors = FALSE, sep = ',')) #load pyscenic output

auc_sub <- auc_mtx[c('CEBPB_', 'NF1_', 'ATF4_', 'FOXJ2_', 'EGR1_',
                     'HLF_', 'CEBPD_', 'NR1H4_', 'CEBPA_',
                     'STAT3_', 'NF1_', 'ELF3_', 'STAT6_', 'NR2F2_',
                     'JUN_', 'FOS_', 'ETS1_', 'PAX5_'),]

auc_sub <- t(auc_sub)

colnames(auc_sub) <- c('CEBPB_Activity', 'NF1_Activity', 'ATF4_Activity', 'FOXJ2_Activity', 'EGR1_Activity',
                       'HLF_Activity', 'CEBPD_Activity', 'NR1H4_Activity', 'CEBPA_Activity',
                       'STAT3_Activity', 'NF1_Activity', 'ELF3_Activity', 'STAT6_Activity', 'NR2F2_Activity',
                       'JUN_Activity', 'FOS_Activity', 'ETS1_Activity', 'PAX5_Activity')

st_hcc06@meta.data <- merge(st_hcc06@meta.data, auc_sub,
                            by = 'row.names', all = TRUE)

rownames(st_hcc06@meta.data) <- st_hcc06@meta.data$Row.names

SpatialFeaturePlot(object = st_hcc06, 
                   features = c('CEBPB_Activity', 'NF1_Activity', 'ATF4_Activity', 'EGR1_Activity',
                                'HLF_Activity', 'CEBPD_Activity', 'NR1H4_Activity', 'CEBPA_Activity',
                                'STAT3_Activity', 'NF1_Activity', 'ELF3_Activity', 'STAT6_Activity',
                                'JUN_Activity', 'FOS_Activity', 'ETS1_Activity', 'PAX5_Activity', 'XBP1_Activity') , 
                   alpha = c(0.1, 1), ncol=4) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")

#-- Module Scores ---------------------------------------------
CAF_list <- list(c('COL6A3', 'COL3A1', 'COL1A1', 'COL1A2' ,'VIM', 'DCN'))
st_hcc06 <- AddModuleScore(object = st_hcc06, features = CAF_list, name = "CAF_score")
SpatialFeaturePlot(object = st_hcc06, features = "CAF_score1", alpha = c(0.1, 1))

BCellList <- list(c('CD19', 'CD79A', 'CD40', 'CR2', 'FCER2A',
                    'CD22', 'TNFRSF13C', 'IL7R', 'BLNK', 'CIITA',
                    'IGHM', 'IGHD', 'IGHG3','IGHG1','IGKC','IGLC1','IGLC2'))
st_hcc06 <- AddModuleScore(object = st_hcc06, features = BCellList, name = "Bcell_score")
SpatialFeaturePlot(object = st_hcc06, features = "Bcell_score1", alpha = c(0.1, 1))


plasemaCellList <- list(c('PRDM1', 'CXCR4','CD27', 'SDC1', 
                          'IGHM', 'IGHD', 'IGHG3','IGHG1',
                          'IGKC','IGLC1','IGLC2', 'IGHA1'))
st_hcc06 <- AddModuleScore(object = st_hcc06, features = plasemaCellList, name = "Plasma_score")
SpatialFeaturePlot(object = st_hcc06, features = "Plasma_score1", alpha = c(0.1, 1))


angioList <- list(c('PECAM1', 'VWF', 'ENG', 'CDH5', 'CD34'))
st_hcc06 <- AddModuleScore(object = st_hcc06, features = angioList, name = "VasScore")
SpatialFeaturePlot(object = st_hcc06, features = "VasScore1", alpha = c(0.1, 1))


#-- Spatial Plots ---------------------------------------------
SpatialFeaturePlot(object = st_hcc06, features = c('ANG',  'FGA', 'HP', 'C3',
                                                   'APOA1', 'AFP', 'VEGFA',
                                                   'ALDH2', 'TGFB1', 'IL6R'), 
                   alpha = c(0.1, 1), ncol = 4) & scale_fill_continuous(low = "white", high = "red")



SpatialFeaturePlot(object = st_hcc06, 
                   features = c('COL1A1', 'COL1A2', 'COL6A3', 'COL6A1', 'COL3A1','VIM', 'THY1', 'DCN') , 
                   alpha = c(0.1, 1), ncol=3) & scale_fill_continuous(low = "white", high = "red")

## subsetting
st_hcc04_sub   <- subset(x = st_hcc06,  idents = c(4, 8))
st_hcc04_tumor <- subset(x = st_hcc06,  idents = c(0, 1, 2))


st_hcc06 <- subset(x = st_hcc06, subset = FOS_Activity < 0.12)
st_hcc06 <- subset(x = st_hcc06, subset = EGR1_Activity < 0.16)

SpatialFeaturePlot(object = st_hcc06, 
                   features = c('NR2F2_Activity', 'FOS_Activity', 
                                'EGR1_Activity') , 
                   alpha = c(0.1, 1), ncol=2) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


SpatialFeaturePlot(object = st_hcc06, 
                   features = c('HLF_Activity', 'CEBPA_Activity', 'NR1H4_Activity') , 
                   alpha = c(0.1, 1), ncol=3) & theme(legend.position='right') &
  scale_fill_continuous(low = "white", high = "green")


VlnPlot(object = st_hcc06_tumor, features = c('CEBPD_Activity', 'NR1H4_Activity', 
                                              'CEBPA_Activity', 'ATF4_Activity',
                                              'STAT3_Activity'),
        split.by = 'seurat_clusters', ncol=2, pt.size=0) & 
  stat_compare_means(size=5)  & 
  theme(text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')

scatterPlot(st_hcc06_sub, 'STAT3', 'IL6R', 'IL6R') 



SpatialFeaturePlot(object = st_hcc06, features = c('EPCAM', 'THY1', 'KRT19', 'PROM1',
                                                   'ALDH1A1', 'CD24', 'ANPEP','CD44',
                                                   'ICAM1', 'CD47', 'SOX9', 'ALDH1A1',
                                                   'ALDH2', 'ABCG2'), 
                   alpha = c(0.1, 1), ncol = 3) & scale_fill_continuous(low = "white", high = "red")& 
  theme(legend.position='right')


SpatialFeaturePlot(object = st_hcc06, features = c('VEGFA', 'ANG', 'ANGPTL3',
                                                   'HLF_Activity', 'CEBPA_Activity', 'NR1H4_Activity'), 
                   alpha = c(0.1, 1), ncol = 3) & scale_fill_continuous(low = "white", high = "red") & 
  theme(legend.position='right')


VlnPlot(object = st_hcc06, features = c('VasScore1'),
        split.by = 'seurat_clusters') & 
  stat_compare_means(size=6)  & 
  ggtitle('Angiogenic Factor Score') &
  theme(text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')


########################## Angiokine ###########################

angiokine <- c('ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
               'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
               'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
               'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')

angiokine_573 <- c( 'ANGPT1', 'ANGPT2', 'AXL', 'FGF2', 'KIT', 'CA9', 'EPO',
                    'FGF19', 'FGF21', 'FGF23', 'GAS6', 'HGF', 'MET', 'PDGFB',
                    'PIGF', 'KRAS', 'HRAS', 'NRAS', 'RAS', 'VEGFA', 'VEGFB',
                    'VDGFC', 'VEGFD', 'FLT1', 'KDR', 'FLT4')



VlnPlot(object = st_hcc06, features = angiokine,
        split.by = 'seurat_clusters', ncol=6, pt.size=0) & 
  stat_compare_means(size=4)  & 
  theme(text = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) & 
  theme(legend.position='right')


#####################################          GESA              ###########################

# Find highly expressed genes in each cluster
hcc06_markers <- FindAllMarkers(st_hcc06, 
                                test.use = "negbinom", 
                                only.pos = T, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25)
markers_df <- as.data.frame(hcc06_markers %>% 
                              group_by(cluster) %>% 
                              slice_max(n = 5, order_by = avg_log2FC))

library(clusterProfiler)
library(org.Hs.eg.db)
library(EnrichmentBrowser)
library(AnnotationDbi)
# Extract the foldchanges for GSEA
res <- hcc06_markers[hcc06_markers$cluster == 2,]
foldchanges <- res$avg_log2FC
names(foldchanges) <- res$gene
#names(foldchanges) <- rownames(res)
foldchanges <- sort(foldchanges, decreasing = TRUE)
# use mapIds method to obtain Entrez IDs
names(foldchanges) <- mapIds(org.Hs.eg.db, names(foldchanges), 'ENTREZID', 'SYMBOL')
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    keyType = "kegg",
                    nPerm = 1000, # default number permutations
                    minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

