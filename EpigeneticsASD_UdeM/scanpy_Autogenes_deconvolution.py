
pip install -U matplotlib

pip install autogenes

import numpy as np
import scanpy as sc
import scipy as sci
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import autogenes as ag

from sklearn.svm import NuSVR
import pickle

#read single-cell data
file = './expreMatrixEnum.csv.gz'
adata = sc.read(file, cache=False).transpose()
adata.var_names

#normalizing and selecting 4000 hihgly variable genes for optimization
#we use log normalized data for selecting hihgly variable genes and visualization
adata_norm = sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4, copy=True) 
adata_log = sc.pp.log1p(adata_norm, copy=True) 
sc.pp.highly_variable_genes(adata_log, flavor='cell_ranger', n_top_genes=4000)
adata_proc = adata_norm[:, adata_log.var[adata_log.var['highly_variable']==True].index]
adata_proc

sc.pp.pca(adata_log, n_comps=30, use_highly_variable=True, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_log, log=True)

adata_log.obs['cells'] = [x.split('.', 1)[0] for x in adata_log.obs_names]

adata_log.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat

sc.pl.pca_scatter(adata_log, color='cells', save='PCA_forebrain_4000top_genes.png')
#plt.savefig(sc.pl.pca_scatter(adata_log, color='cells'), 'PCA_brain_4000top_genes.png')

np.unique([x.split('.', 1)[0] for x in adata_log.obs_names])

#filter cells in normalized data
adata_proc = adata_proc[adata_log.obs_names]

#calculating the centroids of cell types
#change the vector with the corresponding cell names 
clusters = np.array(["Excitatory neurons possibly midbrain","Excitatory neurons midbrain","GABAergic forebrain","Forebrain early neuroblast possibly GABAergic","Cortical Interneurons","Forebrain neural progenitor EMX1","Neuroblast motorneuron/GABAergic?","Glioblast","Radial Glia cycling","Radial glia/Glioblast/Forebrain progenitor EMX1","Cortical Pyramidal","Excitatory neurons cortex","Radial Glia potential glioblast","Radial Glia VLMC primed?","OPCs","Radial Glia","Microglia","Glioblast/Pre-OPC","U","Striatum/Cortical neurons","Radial Glia/Glioblast","Forebrain inhibitory neuroblast","Hindbrain neuroblast","Endothelial","VLMCs","Inhibitory neurons Midbrain, possibly GABAergic","Mid/Hindbrain neuroblast","GABAergic or interneuron neuroblast probably midbrain","Midbrain inhibitory neuroblast","Midbrain/Hindbrain inhibitory neuroblast"])
sc_mean = pd.DataFrame(index=adata_proc.var_names,columns=clusters)
for cluster in clusters:
    cells = [x for x in adata_proc.obs_names if x.startswith(cluster)]
    sc_part = adata_proc[cells,:].X.T
    sc_mean[cluster] = pd.DataFrame(np.mean(sc_part,axis=1),index=adata_proc.var_names)
    
centroids_sc_hv = sc_mean
centroids_sc_hv.shape

ag.init(centroids_sc_hv.T)
ag.optimize(ngen=5000,seed=0,nfeatures=400,mode='fixed',offspring_size=100,verbose=False)
ag.plot(weights=(-1,0))

index = ag.select(index=0)

#filter marker genes in the bulk samples
centroids_sc_pareto = centroids_sc_hv[index]

#Correlation matrix
corr = pd.DataFrame(data = np.corrcoef(centroids_sc_pareto.T), columns = centroids_sc_pareto.columns, index = centroids_sc_pareto.columns)
mask = np.zeros_like(corr)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    sns_plot =sns.clustermap(np.abs(corr),cmap=sns.color_palette("GnBu", 100), robust=True)



#marker genes
import seaborn as sns
subTypes = pd.DataFrame
subTypes = centroids_sc_pareto.columns
type_pal = sns.husl_palette(centroids_sc_pareto.columns.size, s=0.7)
lut = dict(zip(centroids_sc_pareto.columns.unique(), type_pal))
row_colors = subTypes.map(lut)
sns_plot = sns.clustermap(centroids_sc_pareto.T, cmap="mako", robust=True)

centroids_sc_pareto

"""# Regression"""

data_bulk_raw = pd.read_csv('expreMatrixEnum.csv',delimiter='\t', index_col=0)


coef_nusvr = ag.deconvolve(data_bulk_raw.T, model='nusvr')
coef_nnls = ag.deconvolve(data_bulk_raw.T, model='nnls')

def normalize_proportions(data,copy):
    if copy==True:
        data_copy = data.copy()
    else:
        data_copy = data
    data_copy[data_copy < 0] = 0
    for raw in data_copy.index:
        sum = data_copy.loc[raw].sum()
        data_copy.loc[raw] = np.divide(data_copy.loc[raw],sum)
    return data_copy

proportions_NuSVR = normalize_proportions(pd.DataFrame(data=coef_nusvr,columns=clusters,index=data_bulk_raw.columns), copy = False)
proportions_nnls = normalize_proportions(pd.DataFrame(data=coef_nnls,columns=clusters,index=data_bulk_raw.columns), copy = False)

proportions_nnls
#centroids_sc_pareto.index

proportions_nnls.to_csv("proportions_nnls_Interneurons_E18.5_6000genes_200markers.tsv")
proportions_NuSVR.to_csv("proportions_NuSVR_Interneurons_E18.5_6000genes_200markers.tsv")
