#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 00:23:30 2023

@author: chenxuanwei
"""

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt

import time
import pickle as save
from memory_profiler import profile
import os
#data input
fn=1
os.chdir('./'+str(fn))
name=["brainA"]
k=0
dfindex=str(name[k])+'_counts.csv'
df=pd.read_csv(dfindex, header=0, index_col=0)
df = df.T[df.sum(0) >= 3].T 
zero_counts = (df == 0).sum().sum()
zero_ratio = zero_counts/(df.shape[0]*df.shape[1])
# Align count matrix with metadata table
df['total_counts'] = df.sum(axis=1)

for i in range(df.shape[0]):
    str_1 = df.index[i][0:]
    str_spit = str_1.split(sep='x')
    df.at[df.index[i],'x'] = float(str_spit[0])
    df.at[df.index[i],'y'] = float(str_spit[1])
    
    
sample_info = df[['x','y','total_counts']]
df.drop(['x','y','total_counts'],axis=1,inplace=True)
#plt.scatter(sample_info['x'], sample_info['y'], c=sample_info['total_counts'])
#plt.axis('equal')
    
#sample_info = sample_info.query('total_counts > 10')  
df = df.loc[sample_info.index]


#SOMDE
import NaiveDE
import SpatialDE
corinfo = sample_info
X = corinfo[['x','y']].values.astype(np.float32)
df = df.loc[sample_info.index]

from somde import SomNode
@profile
def count():
    som = SomNode(X,10)

    ndf,ninfo = som.mtx(df.T)

    dfm = NaiveDE.stabilize(ndf)

    dfm = dfm.fillna(0)
    som.nres = NaiveDE.regress_out(ninfo, dfm, 'np.log(total_count)').T

    result,SVnum = som.run()
    
    filename = str(name[k])+"_result_SOMDE.data"
    resultname = str(name[k])+"_result_SOMDE.csv"
    result.to_csv(resultname)

if __name__ == '__main__':
    count()


    


