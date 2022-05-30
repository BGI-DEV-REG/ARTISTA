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
data = sc.read('Adult_telencephalon_rep2_DP8400015234BLA3_1_SCT_Removed_0305.h5ad')

### Spatial clustering
sc.pp.pca(data,n_comps=20)
sc.pp.neighbors(data, n_neighbors=30)

sq.gr.spatial_neighbors(data, n_neighs = 6)
data.obsp['connectivities'].data
conn = data.obsp['connectivities'].copy()
conn.data[conn.data > 0] = 1
data.obsp['spatial_distances'].data
adj = conn + data.obsp['spatial_connectivities']
adj.data[adj.data  > 0] = 1

sc.tl.umap(data)

### scanpy plot 
sc.tl.leiden(data, adjacency=adj, key_added='spatial_leiden',resolution=0.5)
sc.pl.spatial(data, spot_size = 30, basis = 'spatial', color=['spatial_leiden'])

### save scanpy plot
outpath = 'your path'
idnames = 'Adult_telencephalon_rep2_DP8400015234BLA3_1'
plt.savefig(outpath+'/'+idnames+'.0.5.pdf')

### Color
colorlist=pd.read_table('Adult_telencephalon_rep2_DP8400015234BLA3_1_Brain_Color.txt',sep='\t',header=None)
colorlist.columns = ['order','Color']
colorlist['Color'] =colorlist['Color'].apply(str)
colorlist['order'] =colorlist['order'].apply(str)

### Spatial plot
featureplot_slices_discrete(obj = data,
    feature = 'spatial_leiden',
    fname = os.path.join(outpath,f'{idnames}.0.5.spatial.pdf'),
    show = False,
    scale = True,
    legend_size = 6,
    order = colorlist['order'].tolist(),
    colors = colorlist['Color'].tolist(),
    slices = None,
    angle_dict = data.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    compress_factor=False,
    #dpi = 300,
    raw=False)

