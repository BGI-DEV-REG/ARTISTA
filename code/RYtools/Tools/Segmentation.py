# -*- coding: utf-8 -*-
"""
RY Segmentation

sample:
seg_obj = Segmentation()
seg_obj.load(img_path, mRNA_path)  # ssDNA image and mRNA image
seg_obj.pre_process(threshold=threshold)
seg_obj.watershed(block_size=block_size,
                  offset=offset,
                  min_distance=min_distance)
seg_obj.save_scGEM(save_path, name)
"""


import matplotlib.pyplot as plt
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
import cv2
from scipy import ndimage as ndi
from skimage import (
    color, feature, filters, measure, segmentation, io
)
import numpy as np
import pandas as pd
import os


class Segmentation:

    def __init__(self):
        self.img_path = None
        self.mRNA_path = None

        self.raw_img = None
        self.img = None
        self.mask = None
        self.label = None

    def load(self, img_path, mRNA_path, signal_pbar=None):
        self.mRNA_path = mRNA_path
        self.img_path = img_path
        self.raw_img = cv2.imread(self.img_path, cv2.IMREAD_GRAYSCALE)

    def pre_process(self,
                    threshold='auto',
                    verbose=True,
                    signal_pbar=None
                    ):
        if threshold == 'auto':
            threshold, _ = cv2.threshold(self.raw_img.copy(), 0, 255, cv2.THRESH_OTSU)
        _, self.img = cv2.threshold(self.raw_img.copy(), threshold, 255, cv2.THRESH_TOZERO)
        if verbose:
            print(f'Used Threshold: {threshold}')
            plt.figure(figsize=(16, 16))
            plt.imshow(self.img, 'gray')

    def watershed(self,
                  block_size=41,
                  offset=0.003,
                  min_distance=15,
                  verbose=True,
                  signal_pbar=None
                  ):
        img = self.img.copy()
        threshold = filters.threshold_local(img, block_size=block_size, offset=offset)
        if verbose:
            plt.figure(figsize=(16, 16))
            plt.imshow(threshold)
            plt.title('Local threshold')
        distance = ndi.distance_transform_edt(img > threshold)
        if verbose:
            plt.figure(figsize=(16, 16))
            plt.imshow(distance)
            plt.title('Distance map')

        local_max_coords = feature.peak_local_max(distance, min_distance=min_distance)
        local_max_mask = np.zeros(distance.shape, dtype=bool)
        local_max_mask[tuple(local_max_coords.T)] = True
        markers = measure.label(local_max_mask)
        if verbose:
            plt.figure(figsize=(16, 16))
            plt.imshow(markers)
            plt.title('Markers')

        self.mask = segmentation.watershed(-distance, markers, mask=img)
        label = color.label2rgb(self.mask, bg_label=0)
        self.label = (label * 255).astype(np.uint8)
        if verbose:
            print("numbers of cells:", self.mask.max())
            plt.figure(figsize=(16, 16))
            plt.imshow(self.label)
            plt.title('Label image')

    # def cellpose(self,
    #              gpu=False,
    #              model_type='cyto',
    #              diameter=None,
    #              flow_threshold=0.4,
    #              verbose=True,
    #              signal_pbar=None
    #              ):
    #     from cellpose import models
    #     img = self.img.copy()
    #     model = models.Cellpose(gpu=gpu, model_type=model_type)
    #     masks, flows, styles, diams = model.eval([img], diameter=diameter, channels=[0, 0],
    #                                              flow_threshold=flow_threshold, do_3D=False)
    #     self.mask = masks[0]
    #     label = color.label2rgb(self.mask, bg_label=0)
    #     self.label = (label * 255).astype(np.uint8)
    #     if verbose:
    #         print("numbers of cells:", self.mask.max())
    #         plt.figure(figsize=(16, 16))
    #         plt.imshow(self.label)

    def save_scGEM(self,
                   save_path,
                   name,
                   verbose=True,
                   signal_pbar=None
                   ):
        data = pd.read_csv(self.mRNA_path, sep='\t',comment="#")

        seg_cell_coor = []
        min_x = data['x'].min()
        min_y = data['y'].min()
        for i in range(self.mask.shape[0]):
            for j in range(self.mask.shape[1]):
                c = self.mask[i, j]
                if c:
                    seg_cell_coor.append([i + min_x, j + min_y, c])
        if signal_pbar:
            signal_pbar.emit(70)
        seg_cell_coor = pd.DataFrame(seg_cell_coor, columns=['x', 'y', 'cell'])
        cell_data = pd.merge(data, seg_cell_coor, how='left', on=['x', 'y'])
        cell_data = cell_data.dropna()
        cell_data['cell'] = cell_data['cell'].astype(int)
        # name = os.path.basename(self.mRNA_path)
        # name = os.path.splitext(name)[0]
        mask_fn = os.path.join(save_path, f'{name}_mask.npy')
        np.save(mask_fn, self.mask)
        gem_fn = os.path.join(save_path, f'{name}_scgem.csv.gz')
        cell_data.to_csv(gem_fn, index=False, sep='\t', compression="gzip")
        # coor_fn = os.path.join(save_path, f'{name}.ssDNA_coor.csv')
        # seg_cell_coor.to_csv(os.path.join(save_path, f'{args.i}.ssDNA_coor.csv'), index=False)
        if verbose:
            print(f'segmented mask save path: {mask_fn}')
            print(f'single-cell GEM save path: {gem_fn}')

    def __repr__(self):
        t = f"ssDNA Image Segmentation Object\n" \
            f"Raw   Image Path: {self.img_path}\n" \
            f"GEM   Data  Path: {self.mRNA_path}"
        return t
