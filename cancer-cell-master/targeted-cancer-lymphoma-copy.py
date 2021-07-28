#!/usr/bin/env python
# coding: utf-8



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns 
import sklearn
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import DistanceMetric
from kneed import KneeLocator as kl
import scanpy as sc
import streamlit as st

print('painpeko')


# sc.settings.verbosity = 2
# sc.logging.print_header()
# sc.settings.set_figure_params(dpi=80, facecolor='lightyellow')

# path = '/Users/devpatelio/Downloads/Coding/Computational_Biology/SingleCell-RNAseq/cancer-cell-master/filtered_feature_bc_matrix'
# results_file = path + 'hodgkins.h5ad'
# adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)


# # In[259]:


# sc.pl.highest_expr_genes(adata, n_top=20)
# print(adata)


# # In[187]:


# #Looking for genes originating from mitochondrial genes
# def percent_mitochondrial(data):
#     data.var['mt'] = data.var_names.str.startswith('MT-')
#     all_keys = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
#     sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=False, log1p=False, inplace=True)
#     diff_percents = set(data.obs.pct_counts_mt.values)
#     print(f"Your dataset consists of cells with the following %s of mitochondrial genes:") 
#     for i in list(diff_percents):
#         print(str(i)+'%')
    
#     return data.obs.pct_counts_mt.value_counts, data.obs.pct_counts_mt, set(data.obs.pct_counts_mt.values)
    
# all_counts_length, all_counts, diff_percents = percent_mitochondrial(adata)

# #Evaluating quality of dataset based on standard measures -> gene counts, total_counts, etc.
# def plot_quality(data, objs, k=3):
#     sc.pl.violin(data, [i for i in objs], multi_panel=True)
#     fig, axes = plt.subplots(1, 3, figsize=(15, 5))
#     sns.distplot(data.obs[objs[0]], ax=axes[0])
#     sns.distplot(data.obs[objs[1]], ax=axes[1])
#     sns.distplot(data.obs[objs[2]], ax=axes[2])

#     sns.jointplot(
#     x=objs[0],
#     y=objs[1],
#     data = data.obs,
#     kind="scatter",
#     s=5
#     )

# plot_quality(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])


# # In[88]:



# def calc_filter_nums(data):
#     sc.pp.filter_cells(data, min_genes=200)
#     sc.pp.filter_genes(data, min_cells=3)

#     sc.pl.scatter(data, x='total_counts', y='pct_counts_mt', color='n_genes_by_counts', title='NoQC_PCT_COUNTS')
#     sc.pl.scatter(data, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', title='NoQC_N_GENES_COUNTS')

#     filter_threshold_gene_counts = int(input('Please enter a gene_count filter slice value:'))
#     filter_threshold_pct_counts = int(input('Please enter a pct_count_mt filter slice value:'))

#     qc_gc = data[data.obs.n_genes_by_counts < filter_threshold_gene_counts, :]
#     sc.pl.scatter(qc_gc, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', title='QC_N_GENES_COUNTS')    
#     try:
#         qc_gc = qc_gc[qc_gc.obs.pct_counts_mt < filter_threshold_pct_counts, :]
#         sc.pl.scatter(qc_gc, x='total_counts', y='pct_counts_mt', color='n_genes_by_counts', title='QC_PCT_COUNTS')
#     except ZeroDivisionError:
#         print('Your total percent of mitochondrial genes are 0%')

#     return qc_gc

# QC_data = adata.copy()
# QC_data = calc_filter_nums(QC_data)


# # In[124]:


# def hvgs (data, target, hyperparams, plot=True):
#     print('*Note, hyperparams should be [min_mean, max_mean, min_disp]')
    
#     def norm_log_data(data, target):
#         sc.pp.normalize_total(data, target_sum=target)
#         sc.pp.log1p(data)

#     norm_log_data(data, target)
#     hvgs = data
#     sc.pp.highly_variable_genes(hvgs, min_mean=hyperparams[0], max_mean=hyperparams[1], min_disp=hyperparams[2])

#     if plot:
#         sc.pl.highly_variable_genes(hvgs)
#     else:
#         pass

#     return hvgs


# norm_data = QC_data.copy()
# HVG = hvgs(norm_data, 1e3, [0.0125, 3, 0.3])


# # In[131]:


# #Remove any effect of the pct_counts after being filtered out and total_counts on data topology (we're looking at highly variable genes so we want to keep the data clean)
# HVGfiltered = HVG[:, HVG.var.highly_variable]
# sc.pp.regress_out(HVGfiltered, ['total_counts', 'pct_counts_mt'])
# sc.pp.scale(HVGfiltered, max_value=10)

# print(len(HVGfiltered))  #HVG Counts
# print(len(adata))  #All Genes Count
# print(f"The total % of all data that is of interest based on its high variability is {round(len(HVG_filtered)/len(adata) * 100, 3)}%") #Calculate total HVG interest


# # In[149]:


# import random

# high_variable_genes_names = HVGfiltered.var_names
# print(f"Total High Variable Genes: {len(high_variable_genes_names)}")

# sc.tl.pca(HVGfiltered, svd_solver='arpack')
# sc.pl.pca(HVGfiltered, color=random.choice(high_variable_genes_names))
# sc.pl.pca_variance_ratio(HVGfiltered, log=True)
# sc.pl.pca_loadings(HVGfiltered)
# HVGfiltered.write(results_file)
# print(HVGfiltered)


# # In[183]:


# def PCA_Elbow_fit(data): #We look each pricniple component based on its proportion of variance, and only include the pcs once a drop in variance is found 
#     model = sklearn.decomposition.PCA()
#     model.fit(data) #calculate PCA on data inputted
#     explained_variance = model.explained_variance_ratio_ #look at the variance ratio or loadings
#     pcs = list(range(1, explained_variance.shape[0]+1)) #renumber all pcs from 1 to all variance_shapes (index 0 provides number)
#     klm = kl(pcs, explained_variance, S=1.0, curve='convex', direction='decreasing') #knee locator finds the drop in variance with a convex curve (can be concave too), direction is based on PC order
#     klm.plot_knee()
#     pcs_used = klm.knee #find the pcs up until the start of the knee (attribute of klm object)
#     pc_list = list(range(1, pcs_used+1))  #lists all numbers of PCs used
#     new_data = sklearn.decomposition.PCA(n_components=pcs_used, svd_solver='arpack').fit_transform(data) #recalculate PCA on data with only the total PCA components found

   
#     return pcs_used, new_data, pcs, explained_variance, pc_list #return number of pcs, PC data embedding with new PCS, number of all initial PCs, variance of previous PCs, list of all new PCs

# label = "HVGfiltered"
# new_frame = pd.DataFrame(HVGfiltered.X, index=HVGfiltered.obs_names, columns=HVGfiltered.var_names) #only want the variable names, objects, and expression values
# pandas_data = new_frame.values #assign values to a new frame


# # In[186]:


# dim, new_pca_data, pc_ax, pc_ay, col_labels = PCA_Elbow_fit(pandas_data) #calculate PCA Elbow fit on formatted data
# print(dim) #17 PCs were used out of the 31

# output_path = path + '_PCA_' + label + str(dim) + '.csv' #output these values to a csv for reference
# PC_frame = pd.DataFrame(new_pca_data, index=new_frame.index.values.tolist(), columns=['PC_' + str(i) for i in col_labels]) #create dataframe with and columns to create dataframe with PC variance ratio with respect to each gene (PC_1 will always have highest variance)
# print(PC_frame.head(2))
# print(PC_frame.shape)


# # In[214]:


# #Computing neighbourhood graphs

# pcaHVG10 = HVGfiltered.copy()

# def tsne_umap_plot(k, var_k, color_n=3, colors=None):
#     if colors == None:
#         colors = random.sample(list(high_variable_genes_names), color_n)
#     sc.pp.neighbors(var_k, n_neighbors=k)   #calculate k-neighbours
    
#     #UMAP
#     sc.tl.umap(var_k)
#     sc.pl.umap(var_k, color=colors, title=f'PCAk={k} UMAP')

#     #TSNE
#     sc.tl.tsne(var_k)
#     sc.pl.tsne(var_k, color=colors, title=f'PCAk={k} tSNE')
    

# tsne_umap_plot(10, pcaHVG10, 4)




# # In[260]:


# #Generating clusters in data
# # fig, axes = plt.subplots(1, 4, figsize=(15, 5))

# sc.tl.leiden(pcaHVG10)
# sc.pl.umap(pcaHVG10, color=['leiden']) #leiden directly uses the k-neighbours plot based on our different PCA embeddings, so the result will look different
# pcaHVG10.write(results_file)


# # In[229]:


# #Labelling clusters with marker genes
# #Utilize 1024 HVGs and their respective names to choose the correct subplot cells based on this cluter method

# sc.tl.rank_genes_groups(pcaHVG10, 'leiden', method='logreg')
# sc.pl.rank_genes_groups(pcaHVG10, n_genes=25, sharey=False)


# # In[254]:


# res = pcaHVG10.uns['rank_genes_groups']

# groups = res['names'].dtype.names
# name_attributes = ['names', 'scores']

# df_markers = pd.DataFrame({g+'-'+key: res[key][g] for g in groups for key in name_attributes})

# #dataframe looks at each cluster and the highest_variable_genes ranked based on its z_score (how similar the hvg to mean expression value of the cluster), higher score indicates this gene's standard deviation is the lowest with respect to the cluster
# marker_names = [i for i in df_markers.iloc[0] if type(i) == str] #get first top row
# print(marker_names)


# # In[261]:


# pcaHVG10.rename_categories('leiden', marker_names)
# sc.pl.umap(pcaHVG10, color='leiden')

# sc.tl.paga(pcaHVG10, groups='leiden')
# sc.pl.paga(pcaHVG10, color=['leiden'] + marker_names)

# sc.pl.stacked_violin(pcaHVG10, marker_names, groupby='leiden')
# sc.pl.dotplot(pcaHVG10, marker_names, groupby='leiden')


# # In[252]:


# #save dataset object
# pcaHVG10.write(results_file, compression='gzip')

