#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt

import time
import pickle as save
from memory_profiler import profile

## Load data after quality control
slideindex='151507'
dfindex=slideindex+'_counts.csv'
df=pd.read_csv(dfindex, header=0, index_col=0)

#df = df.T[df.sum(0) >= 3].T 
#zero_counts = (df == 0).sum().sum()
#zero_ratio = zero_counts/(df.shape[0]*df.shape[1])
# Align count matrix with metadata table
df['total_counts'] = df.sum(axis=1)

for i in range(df.shape[0]):
    str_1 = df.index[i][0:]
    str_spit = str_1.split(sep='x')
    df.at[df.index[i],'x'] = float(str_spit[0])
    df.at[df.index[i],'y'] = float(str_spit[1])
    
sample_info = df[['x','y','total_counts']]
df.drop(['x','y','total_counts'],axis=1,inplace=True)
df = df.loc[sample_info.index]


## SpatialDE
import NaiveDE
import SpatialDE

@profile
def run_wtm_SpatialDE():
    start_time = time.time()
    dfnorm = NaiveDE.stabilize(df.T).T
    res = NaiveDE.regress_out(sample_info, dfnorm.T, 'np.log(total_counts)').T
    res['log_total_count'] = np.log(sample_info['total_counts'])
    X = sample_info[['x', 'y']]
    results = SpatialDE.run(X, res)

    resultname = str(slideindex)+"_result_spatialDE.csv"
    results.to_csv(resultname)
    
    elapsed_time = time.time() - start_time
    et = [elapsed_time]
    filename = str(slideindex)+"_result_spatialDE.data"
    f = open(filename, 'wb')
    save.dump(et, f)

if __name__ == '__main__':
    run_wtm_SpatialDE()

