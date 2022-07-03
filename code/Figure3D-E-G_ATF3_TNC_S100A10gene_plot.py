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


obj = sc.read(r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10/Batch2_Injury_5DPI_rep2_SS200000147BL_D2.h5ad')
out_dir = r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
#max_, min_= 
featureplot_slices_continuous(
    obj = obj,
    cmap='inferno',
    feature = 'AMEX60DDU001003554',
    fname = os.path.join(out_dir, '5DPI_rep2_ATF3_AMEX60DDU001003554.png'),
    show = True,
    scale = True,
    #slices = injs,
    angle_dict = obj.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    min_ = -0.68,
    #max_ = 2,
    #bg_color=(1,1,1,1),
    compress=False,
    raw=False
)

obj = sc.read(r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10/Batch2_Injury_5DPI_rep1_SS200000147BL_D2.h5ad')
out_dir = r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
#max_, min_= 
featureplot_slices_continuous(
    obj = obj,
    cmap='inferno',
    feature = 'AMEX60DDU001003554',
    fname = os.path.join(out_dir, '5DPI_rep1_ATF3_AMEX60DDU001003554.png'),
    show = True,
    scale = True,
    #slices = injs,
    angle_dict = obj.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    min_ = -0.78,
    #max_ = 2,
    #bg_color=(1,1,1,1),
    compress=False,
    raw=False
)


obj = sc.read(r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10/Batch2_Injury_2DPI_rep2_SS200000147BL_D5.h5ad')
obj
out_dir = r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
#max_, min_= 
featureplot_slices_continuous(
    obj = obj,
    cmap='inferno',
    feature = 'AMEX60DDU001003554',
    fname = os.path.join(out_dir, '2DPI_rep2_ATF3_AMEX60DDU001003554.png'),
    show = True,
    scale = True,
    #slices = injs,
    angle_dict = obj.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    min_ = -0.68,
    max_ = 3,
    #bg_color=(1,1,1,1),
    compress=False,
    raw=False
)


obj = sc.read(r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10/Batch2_Injury_2DPI_rep1_SS200000147BL_D5.h5ad')
obj
out_dir = r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
#max_, min_= 
featureplot_slices_continuous(
    obj = obj,
    cmap='inferno',
    feature = 'AMEX60DDU001003554',
    fname = os.path.join(out_dir, '2DPI_rep1_ATF3_AMEX60DDU001003554.png'),
    show = True,
    scale = True,
    #slices = injs,
    angle_dict = obj.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    min_ = -0.68,
    max_ = 2.8,
    #bg_color=(1,1,1,1),
    compress=False,
    raw=False
)

out_dir = r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
#max_, min_= 
featureplot_slices_continuous(
    obj = obj,
    cmap='inferno',
    feature = 'AMEX60DD050822',
    fname = os.path.join(out_dir, '2DPI_rep1_TNC_AMEX60DD050822.png'),
    show = True,
    scale = True,
    #slices = injs,
    angle_dict = obj.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    min_ = -0.4,
    max_ = 3,
    #bg_color=(1,1,1,1),
    compress=False,
    raw=False
)

out_dir = r'../data/Figure3_data/Figure3D-E-G_ATF3_TNC_S100A10'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
#max_, min_= 
featureplot_slices_continuous(
    obj = obj,
    cmap='inferno',
    feature = 'AMEX60DD014975',
    fname = os.path.join(out_dir, '2DPI_rep1_S100A10_AMEX60DD014975.png'),
    show = True,
    scale = True,
    #slices = injs,
    angle_dict = obj.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    min_ = -0.4,
    max_ = 2.5,
    #bg_color=(1,1,1,1),
    compress=False,
    raw=False
)

