import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, cm
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
from Tool.RYSpatial import *
from tqdm import tqdm
import os
import cv2
import numpy as np

l1 = ['Batch1_Injury_15DPI_C_rep1_FP200000266TR_E6','Batch1_Injury_15DPI_rep2_FP200000266TR_E2','Batch1_Injury_15DPI_rep3_FP200000266TR_E3','Batch1_Injury_15DPI_rep4_FP200000266TR_E4','Batch1_Injury_30DPI_rep2_FP200000264BL_A6','Batch1_Injury_60DPI_rep3_FP200000264BL_A6','Batch2_Injury_10DPI_rep1_SS200000147BL_B5','Batch2_Injury_10DPI_rep2_SS200000147BL_B2','Batch2_Injury_10DPI_rep4_SS200000147BL_B3','Batch2_Injury_20DPI_rep1_SS200000147BL_B4','Batch2_Injury_20DPI_rep2_SS200000147BL_B4','Batch2_Injury_20DPI_rep3_SS200000147BL_B5','Batch2_Injury_2DPI_rep1_SS200000147BL_D5','Batch2_Injury_2DPI_rep2_SS200000147BL_D5','Batch2_Injury_2DPI_rep3_SS200000147BL_D4','Batch2_Injury_5DPI_rep1_SS200000147BL_D2','Batch2_Injury_5DPI_rep2_SS200000147BL_D2','Batch2_Injury_5DPI_rep3_SS200000147BL_D3']
l2 = ['1','2','3','4','5','6','8','9','10','11','12','13','14','15','16','17','18','19']

from copy import copy
color_map  = cm.get_cmap('viridis')
color_map  = copy(color_map)
color_map.set_under('lightgray')

out_dir = 'your path'
meta = pd.read_csv('../data/Figure4_data/Figure4J/monocle.anno.csv',index_col=0) ## output from Monocle3
for a,b in zip(l1,l2):
    os.chdir('../data/Figure4_data/Figure4J')
    obj = sc.read(f'{a}.h5ad')
    obj.obs.index = obj.obs.index + f'_{b}'
    obj.obs['pseudotime'] = meta['pseudotime']
    obj.obs['pseudotime'] = obj.obs['pseudotime'].replace(np.nan, -2)
    featureplot_slices_continuous(
        obj=obj,
        feature="pseudotime",
        cmap=color_map,
        fname=os.path.join(out_dir, f'{a}_pseudotime.pdf'),
        show=True,
        scale=True,
        slices=None,
        angle_dict=obj.uns['angle_dict'],
        nrow=1,
        ncol=1,
        max_=30,
        min_=0,
        vmin =0,
        dpi=200,
        compress=True,
        #bg_color=(0,0,0,1),
        raw=False
    )
