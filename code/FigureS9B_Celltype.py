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

os.chdir('../data/FigureS9_data/FigureS9B')
l = [
'Batch1_Injury_15DPI_C_rep1_FP200000266TR_E6',
'Batch1_Injury_15DPI_rep2_FP200000266TR_E2',
'Batch1_Injury_15DPI_rep3_FP200000266TR_E3',
'Batch1_Injury_15DPI_rep4_FP200000266TR_E4']

co = pd.read_csv('Inj_24_Adult_Develop_color_0305.txt',sep='\t')

outpath = 'your path'
for i in l:
    data = sc.read(f'{i}_Annotation_0306.h5ad')
    color_df_sample = co[co['order'].isin(list(set(data.obs['Annotation_0306'])))]
    color_df_sample = color_df_sample.sort_values(by='order')
    featureplot_slices_discrete(obj=data,
                                feature='Annotation_0306',
                                fname=os.path.join(outpath, f'{i}.pdf'),
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
                                raw=False)
