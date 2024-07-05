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
import tracemalloc
import os

for fn in range(1,2):
    #data input
    os.chdir('./'+str(fn))
    name=["brainA"]
    k=0
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    import sepal
    import sepal.datasets as d
    import sepal.models as m
    import sepal.utils as ut
    import sepal.family as family
    import sepal.enrich as fea
    ## sepal
    pth = str(name[k])+"_counts.csv"
    # load in the raw data using the RawData class
    raw_data = d.RawData(pth,)
    raw_data.cnt = ut.filter_genes(raw_data.cnt, min_expr=10, min_occur=5)

    data = m.UnstructuredData(raw_data,eps = 0.1)
    times = m.propagate(data,normalize = True,scale =True)

    resultname = str(name[k])+"_result_sepal.tsv"
    times.to_csv(resultname)
    os.chdir('../')


