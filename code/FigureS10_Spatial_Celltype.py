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
from Tools.Spatial import *

inh5ad = sys.argv[1]
outpath = sys.argv[2]
idnames = sys.argv[3]

os.chdir('../data/Figure3_data/Figure3B')
co = pd.read_csv('Inj_19_color_0305.red.csv')
co['order'] = co['order'].astype(str)
data = sc.read('Batch1_Injury_15DPI_rep2_FP200000266TR_E2_Annotation_0306.h5ad') ## the input data are Control,2DPI-1,5DPI-1,10DPI-1,15DPI-3,30DPI,60DPI
#l = ['REAEGC','RIPC1','IMN']
l = ['NPTXEX','CMPN','WNTEGC','SFRPEGC','REAEGC','MCG','IMN']
for i in l:
    color_df_sample = co[co['order'].isin([i])]
    featureplot_slices_discrete(obj=data,
                                feature='Annotation_0306',
                                fname=os.path.join(outpath, f'{idnames}_{i}_celltype.pdf'),
                                show=False,
                                scale=True,
                                legend_size=6,
                                order=color_df_sample['order'].tolist(),
                                colors=color_df_sample['Color'].tolist(),
                                slices=None,
                                angle_dict=data.uns['angle_dict'],
                                nrow=1,
                                ncol=1,
                                compress_factor=False,
                                dpi = 200,
                                raw=False)

