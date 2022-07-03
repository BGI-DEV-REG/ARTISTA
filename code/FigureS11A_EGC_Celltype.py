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
sys.path.append('../')
from Tools.Spatial import *

os.chdir('../data/FigureS11_data/FigureS11A')
co = pd.read_csv('Inj_24_Adult_Develop_color_0305.txt',sep='\t')
ctype = ['REAEGC','SFRPEGC','WNTEGC']
color_df_sample = co[co['order'].isin(ctype)]

data = sc.read('Batch2_Injury_2DPI_rep1_SS200000147BL_D5_Annotation_0306.h5ad')
outpath = 'your path'
idnames = 'Batch2_Injury_2DPI_rep1_SS200000147BL_D5'
featureplot_slices_discrete(obj = data,
    feature = 'Annotation_0306',
    fname = os.path.join(outpath, f'{idnames}_EGC.pdf'),
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
    dpi = 300,
    raw=False)
