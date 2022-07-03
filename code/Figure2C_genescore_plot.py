import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append(r'/mnt/d/02.project/06RYfish/02brain_dev/00.python_spifig/script/')
from RYSpatial import *
from tqdm import tqdm
import os
import cv2



obj = sc.read(r'/mnt/d/02.project/06RYfish/02brain_dev/00.python_spifig/data/inte_all_v3_ref_rpca_SCT_0827.h5ad')

out_dir = r'/mnt/d/02.project/06RYfish/02brain_dev/00.python_spifig/data/Figure_05292022'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
    
dev = [
 'Stage44_telencephalon_rep2_FP200000239BL_E4',
 'Stage54_telencephalon_rep2_DP8400015649BRD6_2',
 'Stage57_telencephalon_rep2_DP8400015649BRD5_1',
 'Injury_control_FP200000239BL_E3',
 'Adult_telencephalon_rep2_DP8400015234BLA3_1',
 'Meta_telencephalon_rep1_DP8400015234BLB2_1',
]

nes_df = pd.read_csv('../data/Figure2_data/Figure2C/modelscore_dev_0817_1.xls', sep='\t')
nes_df

obj.obs[['Cell_Cycle1', 'NSCgene1','Translation1']] = nes_df[['Cell_Cycle1', 'NSCgene1','Translation1']]
out_dir = '../data/Figure2_data/Figure2C'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
featureplot_slices_continuous(
    obj = obj,
    feature = 'NSCgene1',
    fname = os.path.join(out_dir, 'fig2C_NSCgene1_0822.png'),
    show = True,
    scale = True,
    slices = dev,
    angle_dict = obj.uns['angle_dict'],
    max_ = 1,
    min_ = -0.2,
    nrow = 2,
    ncol = 3,
    dpi=600,
    bg_color=(1,1,1,1),
    compress=False,
    raw=False
)


featureplot_slices_continuous(
    obj = obj,
    feature = 'Translation1',
    fname = os.path.join(out_dir, 'fig2C_Translation1_0822.png'),
    show = True,
    scale = True,
    slices = dev,
    angle_dict = obj.uns['angle_dict'],
    #max_ = 1,
    #min_ = 0,
    nrow = 2,
    ncol = 3,
    dpi=100,
    bg_color=(1,1,1,1),
    compress=False,
    raw=False
)


featureplot_slices_continuous(
    obj = obj,
    feature = 'Cell_Cycle1',
    fname = os.path.join(out_dir, 'fig2C_Cell_Cycle1_0822.png'),
    show = True,
    scale = True,
    slices = dev,
    angle_dict = obj.uns['angle_dict'],
    #max_ = 1,
    #min_ = 0,
    nrow = 2,
    ncol = 3,
    dpi=100,
    bg_color=(1,1,1,1),
    compress=False,
    raw=False
)


