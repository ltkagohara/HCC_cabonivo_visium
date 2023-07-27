```python
RUNNING stLearn ON VISIUM CABO/NIVO DATASET
```


```python
# importing modules
import stlearn as st
import numpy as np
import scanpy as sc
import pandas as pd

st.settings.set_figure_params(dpi=120)
```


```python
# reading in the data
data_hcc2R=st.Read10X("/hcc2R/")
data_hcc2R.var_names_make_unique()
data_hcc6NR=st.Read10X("/hcc6NR/")
data_hcc6NR.var_names_make_unique()
data_hcc7NR=st.Read10X("/hcc7NR/")
data_hcc7NR.var_names_make_unique()
data_hcc3R=st.Read10X("/hcc3R/")
data_hcc3R.var_names_make_unique()
data_hcc4R=st.Read10X("/hcc4R/")
data_hcc4R.var_names_make_unique()
data_hcc1R=st.Read10X("/hcc1R/")
data_hcc1R.var_names_make_unique()
data_hcc5NR=st.Read10X("/hcc5NR/")
data_hcc5NR.var_names_make_unique()
```


```python
# pre-processing for gene count table
st.pp.filter_genes(data_hcc2R,min_cells=3)
st.pp.filter_genes(data_hcc6NR,min_cells=3)
st.pp.filter_genes(data_hcc7NR,min_cells=3)
st.pp.filter_genes(data_hcc3R,min_cells=3)
st.pp.filter_genes(data_hcc1R,min_cells=3)
st.pp.filter_genes(data_hcc1R,min_cells=3)
st.pp.filter_genes(data_hcc5NR,min_cells=3)
st.pp.normalize_total(data_hcc2R)
st.pp.normalize_total(data_hcc6NR)
st.pp.normalize_total(data_hcc7NR)
st.pp.normalize_total(data_hcc3R)
st.pp.normalize_total(data_hcc1R)
st.pp.normalize_total(data_hcc1R)
st.pp.normalize_total(data_hcc5NR)
st.pp.log1p(data_hcc2R)
st.pp.log1p(data_hcc6NR)
st.pp.log1p(data_hcc7NR)
st.pp.log1p(data_hcc3R)
st.pp.log1p(data_hcc1R)
st.pp.log1p(data_hcc1R)
st.pp.log1p(data_hcc5NR)
```


```python
# pre-processing for spot image
st.pp.tiling(data_hcc2R)
st.pp.tiling(data_hcc6NR)
st.pp.tiling(data_hcc7NR)
st.pp.tiling(data_hcc3R)
st.pp.tiling(data_hcc1R)
st.pp.tiling(data_hcc1R)
st.pp.tiling(data_hcc5NR)

st.pp.extract_feature(data_hcc2R)
st.pp.extract_feature(data_hcc6NR)
st.pp.extract_feature(data_hcc7NR)
st.pp.extract_feature(data_hcc3R)
st.pp.extract_feature(data_hcc1R)
st.pp.extract_feature(data_hcc1R)
st.pp.extract_feature(data_hcc5NR)

st.em.run_pca(data_hcc2R,n_comps=50)
st.em.run_pca(data_hcc6NR,n_comps=50)
st.em.run_pca(data_hcc7NR,n_comps=50)
st.em.run_pca(data_hcc3R,n_comps=50)
st.em.run_pca(data_hcc1R,n_comps=50)
st.em.run_pca(data_hcc1R,n_comps=50)
st.em.run_pca(data_hcc5NR,n_comps=50)

data_SME_hcc2R = data_hcc2R.copy()
data_SME_hcc6NR = data_hcc6NR.copy()
data_SME_hcc7NR = data_hcc7NR.copy()
data_SME_hcc3R = data_hcc3R.copy()
data_SME_hcc1R = data_hcc1R.copy()
data_SME_hcc1R = data_hcc1R.copy()
data_SME_hcc5NR = data_hcc5NR.copy()
```


```python
# apply stSME to normalise log transformed data
st.spatial.SME.SME_normalize(data_SME_hcc2R, use_data="raw")
st.spatial.SME.SME_normalize(data_SME_hcc6NR, use_data="raw")
st.spatial.SME.SME_normalize(data_SME_hcc7NR, use_data="raw")
st.spatial.SME.SME_normalize(data_SME_hcc3R, use_data="raw")
st.spatial.SME.SME_normalize(data_SME_hcc1R, use_data="raw")
st.spatial.SME.SME_normalize(data_SME_hcc1R, use_data="raw")
st.spatial.SME.SME_normalize(data_SME_hcc5NR, use_data="raw")
data_SME_hcc2R.X = data_SME_hcc2R.obsm['raw_SME_normalized']
data_SME_hcc6NR.X = data_SME_hcc6NR.obsm['raw_SME_normalized']
data_SME_hcc7NR.X = data_SME_hcc7NR.obsm['raw_SME_normalized']
data_SME_hcc3R.X = data_SME_hcc3R.obsm['raw_SME_normalized']
data_SME_hcc1R.X = data_SME_hcc1R.obsm['raw_SME_normalized']
data_SME_hcc1R.X = data_SME_hcc1R.obsm['raw_SME_normalized']
data_SME_hcc5NR.X = data_SME_hcc5NR.obsm['raw_SME_normalized']
st.pp.scale(data_SME_hcc2R)
st.pp.scale(data_SME_hcc6NR)
st.pp.scale(data_SME_hcc7NR)
st.pp.scale(data_SME_hcc3R)
st.pp.scale(data_SME_hcc1R)
st.pp.scale(data_SME_hcc1R)
st.pp.scale(data_SME_hcc5NR)
st.em.run_pca(data_SME_hcc2R,n_comps=50)
st.em.run_pca(data_SME_hcc6NR,n_comps=50)
st.em.run_pca(data_SME_hcc7NR,n_comps=50)
st.em.run_pca(data_SME_hcc3R,n_comps=50)
st.em.run_pca(data_SME_hcc1R,n_comps=50)
st.em.run_pca(data_SME_hcc1R,n_comps=50)
st.em.run_pca(data_SME_hcc5NR,n_comps=50)
```


```python
# K-means clustering on stSME normalised PCA
# Sample hcc2R has 6 clusters 
st.tl.clustering.kmeans(data_SME_hcc2R,n_clusters=6, use_data="X_pca",
    key_added="X_pca_kmeans")
st.pl.cluster_plot(data_SME_hcc2R, use_label="X_pca_kmeans")
# Sample hcc6NR has 9 clusters
st.tl.clustering.kmeans(data_SME_hcc6NR,n_clusters=9, use_data="X_pca",
    key_added="X_pca_kmeans")
st.pl.cluster_plot(data_SME_hcc6NR, use_label="X_pca_kmeans")
# Sample hcc7NR has 8 clusters
st.tl.clustering.kmeans(data_SME_hcc7NR,n_clusters=8, use_data="X_pca",
    key_added="X_pca_kmeans")
st.pl.cluster_plot(data_SME_hcc7NR, use_label="X_pca_kmeans")
# Sample hcc3R has 6 clusters
st.tl.clustering.kmeans(data_SME_hcc3R,n_clusters=6, use_data="X_pca",
    key_added="X_pca_kmeans")
st.pl.cluster_plot(data_SME_hcc3R, use_label="X_pca_kmeans")
# Sample hcc1R has 10 clusters
st.tl.clustering.kmeans(data_SME_hcc1R,n_clusters=10, use_data="X_pca",
    key_added="X_pca_kmeans")
st.pl.cluster_plot(data_SME_hcc1R, use_label="X_pca_kmeans")
# Sample hcc1R has 9 clusters
st.tl.clustering.kmeans(data_SME_hcc1R,n_clusters=9, use_data="X_pca",
    key_added="X_pca_kmeans")
st.pl.cluster_plot(data_SME_hcc1R, use_label="X_pca_kmeans")
# Sample hcc5NR has 4 clusters
st.tl.clustering.kmeans(data_SME_hcc5NR,n_clusters=4, use_data="X_pca",
    key_added="X_pca_kmeans")
st.pl.cluster_plot(data_SME_hcc5NR, use_label="X_pca_kmeans")

data_SME_hcc2R.uns["lr"] = ['PDCD1_CD274',
                          'PDCD1_CD273',
                          'VSIR_SELPLG', 
                          'CTLA4_CD80', 
                          'CTLA4_CD86',
                          'LAG3_HLA-DRA',
                          'LAG3_FGL1', 
                          'HHLA2_KIR3DL3',
                          'HAVCR2_LGALS9']
data_SME_hcc6NR.uns["lr"] = ['PDCD1_CD274',
                          'PDCD1_CD273',
                          'VSIR_SELPLG', 
                          'CTLA4_CD80', 
                          'CTLA4_CD86',
                          'LAG3_HLA-DRA',
                          'LAG3_FGL1', 
                          'HHLA2_KIR3DL3',
                          'HAVCR2_LGALS9']
data_SME_hcc7NR.uns["lr"] = ['PDCD1_CD274',
                          'PDCD1_CD273',
                          'VSIR_SELPLG', 
                          'CTLA4_CD80', 
                          'CTLA4_CD86',
                          'LAG3_HLA-DRA',
                          'LAG3_FGL1', 
                          'HHLA2_KIR3DL3',
                          'HAVCR2_LGALS9']
data_SME_hcc3R.uns["lr"] = ['PDCD1_CD274',
                          'PDCD1_CD273',
                          'VSIR_SELPLG', 
                          'CTLA4_CD80', 
                          'CTLA4_CD86',
                          'LAG3_HLA-DRA',
                          'LAG3_FGL1', 
                          'HHLA2_KIR3DL3',
                          'HAVCR2_LGALS9']
data_SME_hcc1R.uns["lr"] = ['PDCD1_CD274',
                          'PDCD1_CD273',
                          'VSIR_SELPLG', 
                          'CTLA4_CD80', 
                          'CTLA4_CD86',
                          'LAG3_HLA-DRA',
                          'LAG3_FGL1', 
                          'HHLA2_KIR3DL3',
                          'HAVCR2_LGALS9']
data_SME_hcc1R.uns["lr"] = ['PDCD1_CD274',
                          'PDCD1_CD273',
                          'VSIR_SELPLG', 
                          'CTLA4_CD80', 
                          'CTLA4_CD86',
                          'LAG3_HLA-DRA',
                          'LAG3_FGL1', 
                          'HHLA2_KIR3DL3',
                          'HAVCR2_LGALS9']
data_SME_hcc5NR.uns["lr"] = ['PDCD1_CD274',
                          'PDCD1_CD273',
                          'VSIR_SELPLG', 
                          'CTLA4_CD80', 
                          'CTLA4_CD86',
                          'LAG3_HLA-DRA',
                          'LAG3_FGL1', 
                          'HHLA2_KIR3DL3',
                          'HAVCR2_LGALS9']

st.tl.cci.lr(data_SME_hcc2R)
st.pl.het_plot(data_SME_hcc2R, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc6NR)
st.pl.het_plot(data_SME_hcc6NR, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc7NR)
st.pl.het_plot(data_SME_hcc7NR, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc3R)
st.pl.het_plot(data_SME_hcc3R, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc1R)
st.pl.het_plot(data_SME_hcc1R, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc1R)
st.pl.het_plot(data_SME_hcc1R, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc5NR)
st.pl.het_plot(data_SME_hcc5NR, use_het='cci_lr', image_alpha=0.7)

data_SME_hcc2R.uns["lr"] = ['TNFRSF9_TNFSF9', #4-1BB and 4-1BBL
                          'ICOS_ICOSLG', 
                          'CD28_CD80', 
                          'CD28_CD86',
                          'CD27_CD70',
                          'CD40_CD40LG',
                          'TNFRSF18_TNFSF18', #GITR and GITRL
                          'TNFRSF4_TNFSF4']
data_SME_hcc6NR.uns["lr"] = ['TNFRSF9_TNFSF9', #4-1BB and 4-1BBL
                          'ICOS_ICOSLG', 
                          'CD28_CD80', 
                          'CD28_CD86',
                          'CD27_CD70',
                          'CD40_CD40LG',
                          'TNFRSF18_TNFSF18', #GITR and GITRL
                          'TNFRSF4_TNFSF4']
data_SME_hcc7NR.uns["lr"] = ['TNFRSF9_TNFSF9', #4-1BB and 4-1BBL
                          'ICOS_ICOSLG', 
                          'CD28_CD80', 
                          'CD28_CD86',
                          'CD27_CD70',
                          'CD40_CD40LG',
                          'TNFRSF18_TNFSF18', #GITR and GITRL
                          'TNFRSF4_TNFSF4']
data_SME_hcc3R.uns["lr"] = ['TNFRSF9_TNFSF9', #4-1BB and 4-1BBL
                          'ICOS_ICOSLG', 
                          'CD28_CD80', 
                          'CD28_CD86',
                          'CD27_CD70',
                          'CD40_CD40LG',
                          'TNFRSF18_TNFSF18', #GITR and GITRL
                          'TNFRSF4_TNFSF4']
data_SME_hcc1R.uns["lr"] = ['TNFRSF9_TNFSF9', #4-1BB and 4-1BBL
                          'ICOS_ICOSLG', 
                          'CD28_CD80', 
                          'CD28_CD86',
                          'CD27_CD70',
                          'CD40_CD40LG',
                          'TNFRSF18_TNFSF18', #GITR and GITRL
                          'TNFRSF4_TNFSF4']
data_SME_hcc1R.uns["lr"] = ['TNFRSF9_TNFSF9', #4-1BB and 4-1BBL
                          'ICOS_ICOSLG', 
                          'CD28_CD80', 
                          'CD28_CD86',
                          'CD27_CD70',
                          'CD40_CD40LG',
                          'TNFRSF18_TNFSF18', #GITR and GITRL
                          'TNFRSF4_TNFSF4']
data_SME_hcc5NR.uns["lr"] = ['TNFRSF9_TNFSF9', #4-1BB and 4-1BBL
                          'ICOS_ICOSLG', 
                          'CD28_CD80', 
                          'CD28_CD86',
                          'CD27_CD70',
                          'CD40_CD40LG',
                          'TNFRSF18_TNFSF18', #GITR and GITRL
                          'TNFRSF4_TNFSF4']

st.tl.cci.lr(data_SME_hcc2R)
st.pl.het_plot(data_SME_hcc2R, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc6NR)
st.pl.het_plot(data_SME_hcc6NR, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc7NR)
st.pl.het_plot(data_SME_hcc7NR, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc3R)
st.pl.het_plot(data_SME_hcc3R, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc1R)
st.pl.het_plot(data_SME_hcc1R, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc1R)
st.pl.het_plot(data_SME_hcc1R, use_het='cci_lr', image_alpha=0.7)
st.tl.cci.lr(data_SME_hcc5NR)
st.pl.het_plot(data_SME_hcc5NR, use_het='cci_lr', image_alpha=0.7)
```


```python
# louvain clustering on stSME normalised data
st.pp.neighbors(data_SME,n_neighbors=17,use_rep='X_pca')
st.tl.clustering.louvain(data_SME, resolution=1.19)
st.pl.cluster_plot(data_SME,use_label="louvain",show_plot=True)

data_SME_hcc2R.obs['X_pca_kmeans']
data_SME_hcc6NR.obs['X_pca_kmeans']
data_SME_hcc3R.obs['X_pca_kmeans']
data_SME_hcc7NR.obs['X_pca_kmeans']
data_SME_hcc1R.obs['X_pca_kmeans']
data_SME_hcc4R.obs['X_pca_kmeans']
data_SME_hcc5NR.obs['X_pca_kmeans']

data_SME_hcc2R.write('/data_SME_hcc2R.h5ad')
data_SME_hcc6NR.write('/data_SME_hcc6NR.h5ad')
data_SME_hcc3R.write('/data_SME_hcc3R.h5ad')
data_SME_hcc7NR.write('/data_SME_hcc7NR.h5ad')
data_SME_hcc1R.write('/data_SME_hcc1R.h5ad')
data_SME_hcc4R.write('/data_SME_hcc4R.h5ad')
data_SME_hcc5NR.write('/data_SME_hcc5NR.h5ad')

df = pd.DataFrame(data_SME_hcc2R.obs['X_pca_kmeans'])
df = pd.DataFrame(data_SME_hcc6NR.obs['X_pca_kmeans'])
df = pd.DataFrame(data_SME_hcc3R.obs['X_pca_kmeans'])
df = pd.DataFrame(data_SME_hcc7NR.obs['X_pca_kmeans'])
df = pd.DataFrame(data_SME_hcc1R.obs['X_pca_kmeans'])
df = pd.DataFrame(data_SME_hcc4R.obs['X_pca_kmeans'])
df = pd.DataFrame(data_SME_hcc5NR.obs['X_pca_kmeans'])

df.to_csv (r'/data_SME_hcc2R_identity.csv', index = True, header=True)
df.to_csv (r'/data_SME_hcc6NR_identity.csv', index = True, header=True)
df.to_csv (r'/data_SME_hcc3R_identity.csv', index = True, header=True)
df.to_csv (r'/data_SME_hcc7NR_identity.csv', index = True, header=True)
df.to_csv (r'/data_SME_hcc1R_identity.csv', index = True, header=True)
df.to_csv (r'/data_SME_hcc4R_identity.csv', index = True, header=True)
df.to_csv (r'/data_SME_hcc5NR_identity.csv', index = True, header=True)
```
