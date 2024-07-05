#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import time
import pickle as save
import tracemalloc
from memory_profiler import profile

## Load data after quality control
slideindex='151507'
ff=slideindex+'_counts.csv'

## scGCO
from scGCO import *
import matplotlib
import matplotlib.pyplot as plt
## matplotlib inline
import warnings
warnings.filterwarnings('ignore')
unary_scale_factor=100
label_cost=10
algorithm='expansion'
locs, data =read_spatial_expression(ff,sep=',',num_exp_genes=0.01, num_exp_spots=0.01, min_expression=1)
data_norm = normalize_count_cellranger(data)
print('{}_processing: {}'.format(dfindex,data_norm.shape))

num_cores = mp.cpu_count()
if num_cores > math.floor(data_norm.shape[1]/2):
    num_cores = int(math.loor(data_norm.shap[1]/2))

@profile
def run_wtm_scGCO():
    start_time = time.time()
    tracemalloc.start()
    exp= data_norm.iloc[:,0]
    cellGraph= create_graph_with_weight(locs, exp)

    fig, ax= plt.subplots(1,1,figsize=(5,5)) #, dpi=300)
    ax.set_aspect('equal')

    exp= data_norm.iloc[:,0].values
    cellGraph = create_graph_with_weight(locs, exp)
    ax.scatter(locs[:,0], locs[:,1], s=1, color='black')
    for i in np.arange(cellGraph.shape[0]):
        x = (locs[int(cellGraph[i,0]), 0], locs[int(cellGraph[i,1]), 0]) 
        y = (locs[int(cellGraph[i,0]), 1], locs[int(cellGraph[i,1]), 1])     
        ax.plot(x, y, color='black', linewidth=0.5)
        
    #plt.title('CellGraph')
    #Text(0.5, 1.0, 'CellGraph')

    ## Gene expression classification via Gaussian mixture modeling
    gmmDict= multiGMM(data_norm)
    
    ## Run the main scGCO function to identify genes with a non-random spatial variability
    result_df= identify_spatial_genes(locs, data_norm, 
                                                   cellGraph,gmmDict)
    elapsed_time = time.time() - start_time
    current,peak=tracemalloc.get_traced_memory()
    tracemalloc.stop()

    filename = '{}_results_scGCO.data'.format(slideindex)
    result_df.to_csv('{}_results_scGCO.csv'.format(slideindex))
    et = [elapsed_time,peak,num_cores]
    f = open(filename, 'wb')
    save.dump(et, f)

if __name__ == '__main__':
    run_wtm_scGCO()
    


    

