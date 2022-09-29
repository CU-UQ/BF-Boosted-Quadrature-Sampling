#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 09:54:12 2022

@author: nokicheng
"""

import matplotlib.pyplot as plt
import scipy.io
import numpy as np

file1 = scipy.io.loadmat('B_cor_no_30_TD.mat')
file2 = scipy.io.loadmat('B_cor_no_20_HC.mat')
file_list = [file1,file2]
key_name = ['corr_unif', 'corr_lev', 'corr_lev_vol']
title_name = [['Uniform Sampling','Leverage Score Sampling','Leveraged Volume Sampling'],
              ['Uniform Sampling','Leverage Score Sampling','Leveraged Volume Sampling']]

fig,axs = plt.subplots(2,3,figsize=(12,7))
for j in range(2):
    file = file_list[j]
    for i in np.arange(3):
        ax = axs[j,i]
        data = np.square(file[key_name[i]])
        ax.scatter(data[:,0],data[:,1],c='tab:blue')
        print('Correlation: ',np.corrcoef(data[:,0],data[:,1]))
        ax.set_title(title_name[j][i],fontsize=14)
        ax.set_xlabel(r'Low-fidelity $\mu^2$')
        ax.set_ylabel(r'High-fidelity $\mu^2$')
plt.tight_layout()
plt.savefig('B_cor.eps', format='eps')

