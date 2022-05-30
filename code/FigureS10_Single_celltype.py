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

os.chdir('../data/Figure3_data/Figure3B')
l = [
'Batch1_Injury_15DPI_rep2_FP200000266TR_E2',
'Batch1_Injury_30DPI_rep2_FP200000264BL_A6',
'Batch1_Injury_60DPI_rep3_FP200000264BL_A6',
'Batch1_Injury_control_FP200000239BL_E3',
'Batch2_Injury_10DPI_rep1_SS200000147BL_B5',
'Batch2_Injury_20DPI_rep1_SS200000147BL_B4',
'Batch2_Injury_2DPI_rep1_SS200000147BL_D5',
'Batch2_Injury_5DPI_rep1_SS200000147BL_D2']

co = pd.read_csv('Inj_24_Adult_Develop_color_0305.txt',sep='\t')

outpath = 'your path'
for i in l:
    data = sc.read(f'{i}_Annotation_0306.h5ad')
    color_df_sample = co[co['order'].isin(list(set(data.obs['Annotation_0306'])))]
    color_df_sample = color_df_sample.sort_values(by='order')
    featureplot_single_discrete(obj=data,
                                feature='Annotation_0306',
                                fname=os.path.join(outpath, f'{i}.pdf'),
                                show=False,
                                scale=True,
                                legend_size=6,
                                order=color_df_sample['order'].tolist(),
                                colors=color_df_sample['Color'].tolist(),
                                slice =None,
                                angle_dict=data.uns['angle_dict'],
                                nrow=5,
                                ncol=4,
                                raw=False)
