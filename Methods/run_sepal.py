#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import pickle as save
import tracemalloc
from memory_profiler import profile

## Load data after quality control
slideindex='151507'
pth=slideindex+'_counts.csv'

## sepal
import sepal
import sepal.datasets as d
import sepal.models as m
import sepal.utils as ut
import sepal.family as family
import sepal.enrich as fea

pth = str(slideindex)+"_counts.csv"
raw_data = d.RawData(pth,)
raw_data.cnt = ut.filter_genes(raw_data.cnt, min_expr=10, min_occur=5)

@profile
def run_wtm_sepal():
    start_time = time.time()
    tracemalloc.start()
    data = m.UnstructuredData(raw_data,eps = 0.1)
    times = m.propagate(data,normalize = True,scale =True)

    resultname = str(slideindex)+"_result_sepal.tsv"
    times.to_csv(resultname)

    elapsed_time = time.time() - start_time
    current,peak=tracemalloc.get_traced_memory()
    tracemalloc.stop()

    filename = '{}_results_sepal.data'.format(slideindex)
    et = [elapsed_time,peak]
    f = open(filename, 'wb')
    save.dump(et, f)
    
if __name__ == '__main__':
    run_wtm_sepal()




