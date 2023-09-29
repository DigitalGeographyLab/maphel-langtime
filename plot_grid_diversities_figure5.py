# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:55:09 2023

@author: TuoVaisanen-e01
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# read data in
shannorm = pd.read_pickle(r'W:\maphel_langtime\pickles\new_old_neighbourhoods_norm_shannon.pkl')
uniqnorm = pd.read_pickle(r'W:\maphel_langtime\pickles\new_old_neighbourhoods_norm_unique_langs.pkl')
shanpopw = pd.read_pickle(r'W:\maphel_langtime\pickles\new_old_neighbourhoods_popweigh_shannon.pkl')
uniqpopw = pd.read_pickle(r'W:\maphel_langtime\pickles\new_old_neighbourhoods_popweigh_unique_langs.pkl')

# drop neighbourhoods where population groups have left
shanpopw = shanpopw[shanpopw['celltype'] != 'Ceased']
uniqpopw = uniqpopw[uniqpopw['celltype'] != 'Ceased']
shannorm = shannorm[shannorm['type'] != 'Ceased grid cells']
uniqnorm = uniqnorm[uniqnorm['type'] != 'Ceased grid cells']

# set theme
sns.set()
sns.set(font_scale=1.2)

# start plotting
fig, ax = plt.subplots(2,2, figsize=(15,14))
ax = ax.flatten()

# set palette
palette = {'Somali-inhabited':'C1','Estonian-inhabited':'C0','Finnish-inhabited':'C3'}

# plot unique languages
g = sns.lineplot(x='year', y='score', hue='type', style='celltype', data=uniqpopw,
                 ci=99, estimator=np.nanmean, n_boot=3000, palette=palette, ax=ax[0])
g.set(xlabel='', ylabel='Unique languages', ylim=(0,40))
g.set_title('a.', loc='left', fontsize=15)
handles = ax[0].get_legend().legendHandles
ax[0].get_legend().remove()
ax[0].legend(handles, ['Residential neighbourhood', 'Somali-inhabited','Estonian-inhabited', 'Finnish-inhabited',
                    '\nTemporal type','Present every year', 'New grid cells'])
        
# plot normalized unique languages
g = sns.lineplot(x='year', y='unique_langs', hue='gridtype', style='type', data=uniqnorm,
                 ci=99, estimator=np.nanmean, n_boot=3000, palette=palette, ax=ax[1],
                 legend=False)
g.set(xlabel='', ylabel='Normalized unique languages')
g.set_title('b.', loc='left', fontsize=15)

# plot shannon entropy
g = sns.lineplot(x='year', y='score', hue='type', style='celltype', data=shanpopw,
                 ci=99, estimator=np.nanmean, n_boot=3000, palette=palette, ax=ax[2],
                 legend=False)
g.set(xlabel='', ylabel='Shannon entropy', ylim=(0,1.9))
g.set_title('c.', loc='left', fontsize=15)

# plot normalized shannon entropy
g = sns.lineplot(x='year', y='shannon', hue='gridtype', style='type', data=shannorm,
                 ci=99, estimator=np.nanmean, n_boot=3000, palette=palette, ax=ax[3],
                 legend=False)
g.set(xlabel='', ylabel='Normalized Shannon entropy')
g.set_title('d.', loc='left', fontsize=15)

# save the figure
plt.savefig(r'W:\maphel_langtime\plots\figure5_grid_new_old_trajectories.pdf', dpi=300, bbox_inches='tight')
