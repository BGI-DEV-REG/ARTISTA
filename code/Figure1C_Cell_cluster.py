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

os.chdir('../data/Figure1_data/Figure1C')
meta = pd.read_csv('Adult_telencephalon_rep2_DP8400015234BLA3_1_Removed.meta.csv',index_col=1)

### remove low quality cells
meta.index = meta['CellID']
me = meta[meta['Annotation']!='Remove']

### load data
data = sc.read('Adult_telencephalon_rep2_DP8400015234BLA3_1_raw.h5ad')

### extract quality cells
data = data[me.index]

### creat spatial plot data
a = np.load('Adult_telencephalon_rep2_DP8400015234BLA3_1_cell_segmentation.npy')

me['id'] = me['CellID'].apply(lambda x:int(x.split('.')[1])).values
elements = list(set(me['id']))
a = np.where(~(np.isin(a, elements)), 0, a)

idnames = 'Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1'
data.uns[str(idnames)]={}
data.uns[str(idnames)]['seg_cell']=a
data.obs['Batch'] = str(idnames)

data.obs['Annotation'] = me['Annotation']
data.uns['angle_dict']={}
data.uns['angle_dict'][str(idnames)]=-185

data.obs['cell_id'] = me['id']

### save data
data.write('Adult_telencephalon_rep2_DP8400015234BLA3_1.h5ad')

### color
co = pd.read_csv('Adult_telencephalon_rep2_DP8400015234BLA3_1_Celltype_Color.csv')
color_df_sample = co[co['order'].isin(list(set(data.obs['Annotation_0305'])))]
color_df_sample = color_df_sample.sort_values(by='order')

output = 'your path'

### Celltype plot
featureplot_single_discrete(
    obj = data,
    feature = 'Annotation_0305',
    fname = os.path.join(f'{outpath},/Adult_telencephalon_rep2_DP8400015234BLA3_1_single_celltype.pdf'),
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
