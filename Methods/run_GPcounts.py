#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 15:56:53 2022

@author: xwchen
"""
import os
import numpy as np
import pandas as pd
import time
import pickle as save
from memory_profiler import profile

slideindex='01A'
dfindex=str(slideindex)+'_counts.csv'
df=pd.read_csv(dfindex, header=0, index_col=0)

for i in range(df.shape[0]):
    str_1 = df.index[i]
    str_spit = str_1.split(sep='x')
    df.at[str_1,'x'] = float(str_spit[0])
    df.at[str_1,'y'] = float(str_spit[1])

X = df[['x','y']]
df.drop(['x','y'],axis=1,inplace=True)
Y= df.T


import gpflow
import tensorflow as tf
from GPcounts.GPcounts_Module import Fit_GPcounts

@profile
def run_wtm_GPcounts():
    start_time = time.time()

    gene_name = []
    gene_name = Y.index
    likelihood = 'Negative_binomial' 
    gp_counts = Fit_GPcounts(X,Y.loc[gene_name],safe_mode=False)
    log_likelihood_ratio = gp_counts.One_sample_test(likelihood)
    result = gp_counts.calculate_FDR(log_likelihood_ratio)

    elapsed_time = time.time() - start_time
    resultname = str(slideindex)+"_result_GPcounts.csv"
    result.to_csv(resultname)

    filename = str(slideindex)+"_result_GPcounts.data"
    et = [elapsed_time]
    f = open(filename, 'wb')
    save.dump(et, f)
    f.close()
    
if __name__ == '__main__':
    run_wtm_GPcounts()

