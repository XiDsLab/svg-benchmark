#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 18:56:58 2022

@author: xwchen
"""

import cv2
import numpy as np
import os
import pandas as pd

def image_feature(Im, Spot_coord, Spot_radius):
    ## Step1: Image RGB channels to LAB channels
    ##-------------------------------------------------------------
    Im_LAB = cv2.cvtColor(Im,cv2.COLOR_BGR2LAB)
    Im_LChannel = Im_LAB[:,:,0]
    ##-------------------------------------------------------------
    ##Step2: Obtain the feature within each spot
    num_spot = Spot_coord.shape[0]
    ##
    Spot_feature = np.zeros(shape=(num_spot,1))
    ##
    for spot_index in range(num_spot):
        cur_spot_coord = np.array(Spot_coord.iloc[spot_index,[2,3]])
        cur_spot_coord[0] = round(cur_spot_coord[0]) - 1
        cur_spot_coord[1] = round(cur_spot_coord[1]) - 1
        cur_spot_coord = cur_spot_coord.astype(int)
        ##
        Im_sel = Im_LChannel[(cur_spot_coord[0] - Spot_radius):(cur_spot_coord[0] + Spot_radius + 1):1,(cur_spot_coord[1] - Spot_radius):(cur_spot_coord[1] + Spot_radius + 1):1]
        Spot_feature[spot_index,:] = Im_sel.mean(axis=(0, 1))
    ##
    Spot_meta = Spot_coord
    Spot_meta["Lchannel"] =Spot_feature
    return (Spot_meta,Im_LChannel)

##HE data
#os.chdir('./spatial analysis/HE')

name=["FFPE_Human_Breast_Cancer","FFPE_Human_Cervical_Cancer",
      "FFPE_Human_Normal_Prostate","FFPE_Human_Intestinal_Cancer",
      "FFPE_Mouse_Brain","FFPE_Mouse_Kidney",
      "V1_Human_Heart","V1_Adult_Mouse_Brain",
      "V1_Mouse_Kidney","V1_Mouse_Brain_Sagittal_Anterior_Section_2",
      "V1_Mouse_Brain_Sagittal_Posterior_Section_2","V1_Mouse_Olfactory_Bulb"]

k=0
os.chdir('./'+name[k])
##Image Loading
img_name = name[k]+"_image.tif"
Im = cv2.imread(img_name)

##Spot_coord loading
Spot_coord = pd.read_csv(name[k]+"_coord.csv"
                         ,sep=",",header=0)
##Like
#                  Barcode     imagerow     imagecol
# 0     AAACAAGTATCTCCCA-1  7408.119239  8455.000020
# 1     AAACACCAATAACTGC-1  8485.991089  2740.000007

##Spot_radius
Spot_radius = 43 ## Spot's real radius in 10X ST data

##Run the function
image_feature_tuple = image_feature(Im, Spot_coord, Spot_radius)

##Reult
#image_feature_tuple[0]
#                  Barcode     imagerow     imagecol    Lchannel
# 0     AAACAAGTATCTCCCA-1  7408.119239  8455.000020  131.326199
# 1     AAACACCAATAACTGC-1  8485.991089  2740.000007  140.915048
# 2     AAACAGAGCGACTCCT-1  3095.631956  7905.000019  151.790065
# 3     AAACAGCTTTCAGAAG-1  6569.218977  2052.000005  127.396354

#image_feature_tuple[1] ## The L channel of image
cv2.imwrite(name[k]+"_L.jpg",image_feature_tuple[1])
os.chdir('../')
image_feature_tuple[0].to_csv(name[k]+"_Lchannel.csv")
