#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import cv2
import anndata as ad

# i indicate method i
slideindex = '151507'
# number of expert region annotations
n_clusters=7
for i in range(11):
    adata = pd.read_csv(str(slideindex)+'_counts.csv', header=0, index_col=0)
    
    ## Load top 2,000 SVGs identified by method i
    top2k = pd.read_csv(str(slideindex)+'_top2k.csv', header=0, index_col=0)
    adata = adata[top2k.iloc[:, i]]
    
    index = adata.index
    adata = adata.astype('float64')
    adata = ad.AnnData(adata)
    spatial=pd.read_csv(str(slideindex)+'_info.csv',index_col=0)
    spatial.index = index
    adata.obs["x_array"]=spatial.iloc[:, 0]
    adata.obs["y_array"]=spatial.iloc[:, 1]
    
    ## Select captured samples
    adata.var_names=[i.upper() for i in list(adata.var_names)]
    adata.var["genename"]=adata.var.index.astype("str")
    adata.var_names_make_unique()
    #adata.write_h5ad("151507_data.h5ad")
    
    ## Read in gene expression and spatial location
    #adata=sc.read_h5ad("151507_data.h5ad")

    ## Set coordinates
    x_array=adata.obs["x_array"].tolist()
    y_array=adata.obs["y_array"].tolist()
    adj=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    #spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
    #spg.prefilter_specialgenes(adata)
    ## Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    p=0.5
    ## Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=200)
    print(l)
    

    
    
    ## Set seed
    r_seed=t_seed=n_seed=100
    ## Seaech for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=1.0, step=0.05, tol=5e-3, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    
    clf=spg.SpaGCN()
    clf.set_l(l)
    ## Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    ## Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    
    ## Do cluster refinement(optional)
    ## shape="hexagon" for Visium data, "square" for ST data.
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    
    ## Save results
    df = pd.DataFrame({'x_array': adata.obs['x_array'], 'y_array': adata.obs['y_array'],
                       'pred':adata.obs["pred"], 'refined_pred':adata.obs["refined_pred"]})
    df.to_csv(str(slideindex)+"_"+str(top2k.iloc[:, i].name)+"_result2k.csv")






