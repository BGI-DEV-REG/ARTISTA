import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import time
import os
os.environ["OPENCV_SHRT_MAX"] = str(pow(2,40))
import cv2
sys.path.append('../')
from Tools.Spatial import *

os.chdir('../data/Figure1_data/Figure1C')

data = sc.read('Adult_telencephalon_rep2_DP8400015234BLA3_1_SCT_Removed_0305.h5ad')
co = pd.read_csv('Inj_24_Adult_Develop_color_0305.txt',sep='\t')
#list(set(data.obs['Annotation']).difference(set(old_co['order'])))
color_df_sample = co[co['order'].isin(list(set(data.obs['Annotation_0305'])))]
color_df_sample = color_df_sample.sort_values(by='order')

featureplot_single_discrete(
    obj = data,
    feature = 'Annotation_0305',
    fname = os.path.join('Adult_telencephalon_rep2_DP8400015234BLA3_1_single_celltype.pdf'),
    show = True,
    scale = True,
    legend_size = 6,
    order = color_df_sample['order'].tolist(),
    colors = color_df_sample['Color'].tolist(),
    slice = 'Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1',
    angle_dict = data.uns['angle_dict'],
    nrow = 4,
    ncol = 4,
    compress=True,
    raw=False
)

