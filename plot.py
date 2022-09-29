#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 10:24:46 2022

@author: nokicheng
"""

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io

file1 = scipy.io.loadmat('B_no_18_TD.mat')
file2 = scipy.io.loadmat('B_no_30_TD.mat')
file3 = scipy.io.loadmat('B_no_11_HC.mat')
file4 = scipy.io.loadmat('B_no_18_HC.mat')

key_name = ['err_unif', 'err_unif_boosted','err_lev_score', 'err_lev_score_boosted', 'err_lev_vol', 'err_lev_vol_boosted' ]
name_list = ['Uniform','Uniform','Leverage Score', 'Leverage Score', 'Leveraged Volume', 'Leveraged Volume']

data1 = np.array([file1[name].squeeze() for name in key_name])
err_data1 = data1.reshape(6000)
err_name = np.repeat(name_list,1000)
err_boost = np.repeat(['Regular','Boosted','Regular','Boosted','Regular','Boosted'],1000)
dic = {'Error':err_data1,'Sampling Method':err_name,'boosted':err_boost}
dataf1 = pd.DataFrame(dic,index = np.arange(6000))

data2 = np.array([file2[name].squeeze() for name in key_name])
err_data2 = data2.reshape(6000)
err_name = np.repeat(name_list,1000)
err_boost = np.repeat(['Regular','Boosted','Regular','Boosted','Regular','Boosted'],1000)
dic = {'Error':err_data2,'Sampling Method':err_name,'boosted':err_boost}
dataf2 = pd.DataFrame(dic,index = np.arange(6000))

data3 = np.array([file3[name].squeeze() for name in key_name])
err_data3 = data3.reshape(6000)
err_name = np.repeat(name_list,1000)
err_boost = np.repeat(['Regular','Boosted','Regular','Boosted','Regular','Boosted'],1000)
dic = {'Error':err_data3,'Sampling Method':err_name,'boosted':err_boost}
dataf3 = pd.DataFrame(dic,index = np.arange(6000))

data4 = np.array([file4[name].squeeze() for name in key_name])
err_data4 = data4.reshape(6000)
err_name = np.repeat(name_list,1000)
err_boost = np.repeat(['Regular','Boosted','Regular','Boosted','Regular','Boosted'],1000)
dic = {'Error':err_data4,'Sampling Method':err_name,'boosted':err_boost}
dataf4 = pd.DataFrame(dic,index = np.arange(6000))

fig, axes = plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
fg = sns.boxplot(y='Error',x='Sampling Method',hue='boosted',data=dataf1,ax=axes[0,0])
fg.axhline(file1['err_QR'][0,0],c='b',alpha=0.5)
fg.axhline(file1['err_exact'][0,0],c='y',alpha=0.5)
fg.set_title('Total Degree; m=18',fontsize=14)
fg.set_yscale('log')
fg.legend(title=None,loc=1)

fg = sns.boxplot(y='Error',x='Sampling Method',hue='boosted',data=dataf2,ax=axes[0,1])
fg.axhline(file2['err_QR'][0,0],c='b',alpha=0.5)
fg.axhline(file2['err_exact'][0,0],c='y',alpha=0.5)
fg.set_title('Total Degree; m=30',fontsize=14)
fg.set_yscale('log')
fg.legend(title=None,loc=1)

fg = sns.boxplot(y='Error',x='Sampling Method',hue='boosted',data=dataf3,ax=axes[1,0])
fg.axhline(file3['err_QR'][0,0],c='b',alpha=0.5)
fg.axhline(file3['err_exact'][0,0],c='y',alpha=0.5)
fg.set_title('Hyperbolic Cross; m=12',fontsize=14)
fg.set_yscale('log')
fg.legend(title=None,loc=1)

fg = sns.boxplot(y='Error',x='Sampling Method',hue='boosted',data=dataf4,ax=axes[1,1])
fg.axhline(file4['err_QR'][0,0],c='b',alpha=0.5)
fg.axhline(file4['err_exact'][0,0],c='y',alpha=0.5)
fg.set_title('Hyperbolic Cross; m=20',fontsize=14)
fg.set_yscale('log')
fg.legend(title=None,loc=1)

fig.savefig('beam.eps',format='eps')


