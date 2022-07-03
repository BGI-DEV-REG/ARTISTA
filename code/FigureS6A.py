import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import pandas as pd
import sys
import time
import os
os.environ["OPENCV_SHRT_MAX"] = str(pow(2,40))
import cv2
import dynamo as dyn
sys.path.append('../Project/Axolotl_Brain_Spatial/42.Registration')
from Tools.Spatial import *
inh5ad = sys.argv[1]
idnames = sys.argv[2]


data = sc.read(inh5ad)
outpath='../FigureS6'

color_df = pd.read_csv('../FigureS6/'+ idnames +'_Color.seq.txt', sep='\t')
color_df = color_df.sort_values(by='order')
color_df['order'] = color_df['order'].astype('str')
ct = list(set(data.obs['spatial_leiden_e30_s8']))
color_df_sample = color_df[color_df['order'].isin(ct)]

featureplot_slices_discrete(obj = data,
    feature = 'spatial_leiden_e30_s8',
    fname = os.path.join(outpath,str(idnames)+'.Region_color.pdf'),
    show = False,
    scale = True,
    legend_size = 6,
    order = color_df_sample['order'].tolist(),
    colors = color_df_sample['Color'].tolist(),
    slices = None,
    angle_dict = data.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    compress_factor=False,
    raw=False)
