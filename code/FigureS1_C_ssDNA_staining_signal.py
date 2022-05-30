import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import KernelDensity
import cv2
from scipy.stats import ks_2samp
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import math


def get_data(data, center, h,H=50*2):
    expr = (data['x']>=center[0]-h/2) & (data['x']<center[0]+h/2) & (data['y']>=center[1]-h/2) & (data['y']<center[1]+h/2)
    data = data[expr]
    data['x'] -= center[0]-H//2
    data['y'] -= center[1]-H//2
    return data


def get_c(npy,h,H=50*2):
    contours = []
    for id_ in np.unique(npy):
        if id_ == -1:
            continue
        arr = np.zeros_like(npy, dtype=np.uint8)
        arr[npy==id_] = 255
        contours_, hierarchy = cv2.findContours(arr,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_TC89_L1)
        c = contours_[0] - np.array([(h-H)//2, (h-H)//2])
        contours.append(c)
    return contours


def kde(df, l=[5,10,15,20], t=99.99,h=50*2, kernel='tophat',w=600j,adj=1, name=''):
    cbs = df.copy()
    for i in l:
        kde = KernelDensity(kernel=kernel, 
                        bandwidth=i).fit(X=cbs[['x', 'y']], 
                                          sample_weight=cbs['MIDCounts'])
        xx, yy = np.mgrid[0:h:w, 0:h:w]
        log_d = kde.score_samples(np.vstack([xx.ravel(), yy.ravel()]).T)    
        m = np.exp(log_d.reshape(xx.shape))
        n = m.copy()
        if t < 100:
            m = np.clip(m, m.min(), np.percentile(m, t))
        m -= m.min()
        m /= m.max()*adj
        m = (m*255).astype(np.uint8)
        mm = cv2.applyColorMap(m, cv2.COLORMAP_JET)
        cv2.imwrite(f'FigureS1_C_{kernel}_mm{i}_{name}.png', mm)
        print(f'FigureS1_C_{kernel}_mm{i}_{name}.png')

obj = sc.read('data/Injury_15DPI_rep3_FP200000266TR_E3/Injury_15DPI_rep3_FP200000266TR_E3.ssDNA.0813.h5ad')

# read bin1 data
df = pd.read_csv("data/Injury_15DPI_rep3_FP200000266TR_E3/Injury_15DPI_rep3_FP200000266TR_E3.Gene_Expression_table.tsv.gz",
                 sep='\t')
m = df['x'].min()
n = df['y'].min()
del df

exon = pd.read_csv('data/Injury_15DPI_rep3_FP200000266TR_E3/EXONIC.txt_Injury_15DPI_rep3_FP200000266TR_E3.pick.CB.txt', 
                   sep='_', header=None, names=['x', 'y'])
exon['x'] -= m
exon['y'] -= n
exon['MIDCounts'] = 1
exon = exon.groupby(['x', 'y']).count().reset_index()
intro = pd.read_csv('data/Injury_15DPI_rep3_FP200000266TR_E3/INTRONIC.txt_Injury_15DPI_rep3_FP200000266TR_E3.pick.CB.txt', 
                   sep='_', header=None, names=['x', 'y'])
intro['x'] -= m
intro['y'] -= n
intro['MIDCounts'] = 1
intro = intro.groupby(['x', 'y']).count().reset_index()

# ssDNA image
ssdna = cv2.imread('data/Injury_15DPI_rep3_FP200000266TR_E3/Injury_15DPI_rep3_FP200000266TR_E3_ssDNA.jpg')

# 2D-KDE
h = 2*80
H = 2*50
for center in [(4136, 3704), (4325,3750)]:
    exon_sub = get_data(exon, center, h)
    intro_sub = get_data(intro, center, h)
    npy = obj.uns['seg_cell'][center[0]-h//2:center[0]+h//2, center[1]-h//2:center[1]+h//2]
    contours = get_c(npy, h)
    kde(intro_sub, name=f'intro_{center[0]}_{center[1]}', l=[5], kernel='gaussian', w=100j,h=H)
    kde(exon_sub, name=f'exon_{center[0]}_{center[1]}', l=[5], kernel='gaussian', w=100j,h=H)
    
    # draw contours
    arr =  np.zeros((H, H, 4), dtype=np.uint8)
    arr = cv2.drawContours(arr,contours,-1,(255,255,255,255), 2)
    cv2.imwrite(f'FigureS1_C_contours_{center[0]}_{center[1]}.png', arr)
    ssdna_sub = ssdna[center[0]-H//2:center[0]+H//2, center[1]-H//2:center[1]+H//2]
    cv2.imwrite(f'FigureS1_C_ssdna_{center[0]}_{center[1]}.png', ssdna_sub)