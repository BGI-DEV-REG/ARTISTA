import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append(r'../script')
from RYSpatial import *
from tqdm import tqdm
import os
import cv2
obj = sc.read(r'../h5ad')
obj
out_dir = r'../FigureS6'
genes = pd.read_csv('../FigureS6/marker.txt', sep='\t')
genes = genes['genename'].tolist()
for gene in genes:
    print (gene)
    featureplot_slices_continuous(
        obj = obj,
        #cmap='inferno',
        feature = gene,
        fname = os.path.join(out_dir, f'Exp_{gene}.png'),
        show = True,
        scale = True,
        #slices = injs,
        angle_dict = obj.uns['angle_dict'],
        nrow = 1,
        ncol = 1,
    #    min_ = 0,
    #    max_ = 0.28,
    #bg_color=(1,1,1,1),
        compress=False,
        raw=False
)
    plt.close('all')
