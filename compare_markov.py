# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 10:40:19 2023

@author: TuoVaisanen-e01
"""
import pandas as pd
from scipy.spatial import distance
from scipy.special import kl_div, rel_entr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# read in population probability matrices
est = pd.read_pickle(r'W:\maphel_langtime\pickles\markov_probs_estpop_shannon_natbre_k5.pkl')
som = pd.read_pickle(r'W:\maphel_langtime\pickles\markov_probs_sompop_shannon_natbre_k5.pkl')
#hma = pd.read_pickle(r'W:\maphel_langtime\pickles\markov_probs_pop_count_shannon_jenks_k5.pkl')
fin = pd.read_pickle(r'W:\maphel_langtime\pickles\markov_probs_finpop_shannon_natbre_k5.pkl')
foreign = pd.read_pickle(r'W:\maphel_langtime\pickles\markov_probs_foreign_pop_shannon_natbre_k5.pkl')

# convert dataframes to numpy arrays
ahma = hma.to_numpy()
aest = est.to_numpy()
asom = som.to_numpy()
afin = fin.to_numpy()
afor = foreign.to_numpy()

# calculate JS distances, values are transposed so it compares correct values
somdif = distance.jensenshannon(ahma.T, asom.T)
estdif = distance.jensenshannon(ahma.T, aest.T)
sedif = distance.jensenshannon(aest.T,asom.T)
somdif2 = distance.jensenshannon(afin.T, asom.T)
estdif2 = distance.jensenshannon(afin.T, aest.T)
natfor = distance.jensenshannon(afin.T, afor.T)
somfor = distance.jensenshannon(asom.T, afor.T)
estfor = distance.jensenshannon(afor.T, aest.T)

# generate array of js distances
js = np.array([somdif, estdif, sedif])
js2 = np.array([somdif2, estdif2, sedif])
js3 = np.array([somfor, estfor, natfor])


# plot the matrices
sns.set()
fig, axes = plt.subplots(2,2, figsize=(14,12))

# set tick labels
ticks = ['Very Low', 'Low', 'Moderate', 'High', 'Very High']

# flatten axes for easier looping
axes = axes.flatten()

# plot finnish speakers
g = sns.heatmap(afin, annot=True, linewidths=.1, ax=axes[0], cbar=True,
                vmin=0, vmax=1, square=True, cmap='OrRd', fmt='.3f',
                xticklabels=ticks, yticklabels=ticks, cbar_kws={'shrink':0.9})
g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
axes[0].set_title('a.', fontsize=18, loc='left')

# plot somali speakers
g = sns.heatmap(asom, annot=True, linewidths=.1, ax=axes[1], cbar=True,
                vmin=0, vmax=1, square=True, cmap='YlGn', fmt='.3f',
                xticklabels=ticks, yticklabels=ticks, cbar_kws={'shrink':0.9})
g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
axes[1].set_title('b.', fontsize=18, loc='left')

# plot estonian speakers
g = sns.heatmap(aest, annot=True, linewidths=.1, ax=axes[2], cbar=True,
                vmin=0, vmax=1, square=True, cmap='YlGnBu', fmt='.3f',
                xticklabels=ticks, yticklabels=ticks, cbar_kws={'shrink':0.9})
g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
axes[2].set_title('c.', fontsize=18, loc='left')

# plot forlang speakers
g = sns.heatmap(afor, annot=True, linewidths=.1, ax=axes[3], cbar=True,
                vmin=0, vmax=1, square=True, cmap='RdPu', fmt='.3f',
                xticklabels=ticks, yticklabels=ticks, cbar_kws={'shrink':0.9})
g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
axes[3].set_title('d.', fontsize=18, loc='left')

plt.savefig(r'W:\maphel_langtime\plots\markov_global_comparison_popgroup_shannon_natbre.pdf', dpi=300,
            bbox_inches='tight')

# get table of similarities
jsdf = pd.DataFrame(js, index=['Somali-HMA', 'Estonian-HMA', 'Somali-Estonian'],
                    columns=ticks)

# get tables of similarities
jsdf2 = pd.DataFrame(js2, index=['Somali-Finnish', 'Estonian-Finnish', 'Somali-Estonian'],
                    columns=ticks)
jsdf3 =  pd.DataFrame(js3, index=['Somali-Foreign', 'Estonian-Foreign','Finnish-Foreign'],
                      columns=ticks)

# plot the distance matrix
fig, ax = plt.subplots(figsize=(7,4))
g = sns.heatmap(jsdf, annot=True, fmt='.3f', square=True, cmap='mako_r',
                vmin=0, vmax=1, cbar=True, cbar_kws={'shrink':0.8}, ax=ax)
plt.savefig(r'W:\maphel_langtime\plots\jensen-shannon_distances_hma.pdf', dpi=300,
            bbox_inches='tight')

# plot the distance matrix
fig, ax = plt.subplots(1, 2, figsize=(15,5))
ax = ax.flatten()
g = sns.heatmap(jsdf2, annot=True, fmt='.3f', square=True, cmap='mako_r',
                vmin=0, vmax=1, cbar=True, cbar_kws={'shrink':0.7}, ax=ax[0])
ax[0].set_title('a.', loc='left', fontsize=16)
g = sns.heatmap(jsdf3, annot=True, fmt='.3f', square=True, cmap='mako_r',
                vmin=0, vmax=1, cbar=True, cbar_kws={'shrink':0.7}, ax=ax[1])
ax[1].set_title('b.', loc='left', fontsize=16)
plt.savefig(r'W:\maphel_langtime\plots\jensen-shannon_distances_pop_group_shannon_naive_quint.pdf', dpi=300,
            bbox_inches='tight')

# plot global markov comparisons for nat and foreigns as well
fig, axes = plt.subplots(2,2, figsize=(15,11))

# flatten axes for easier looping
axes = axes.flatten()

natcolor = sns.cubehelix_palette(start=0.1, rot=.65, dark=0.25, light=.97, as_cmap=True)

# plot native population
g = sns.heatmap(afise, annot=True, linewidths=.1, ax=axes[0], cbar=True,
                vmin=0, vmax=1, square=True, cmap=natcolor, fmt='.3f',
                xticklabels=ticks, yticklabels=ticks, cbar_kws={'shrink':0.95})
g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
axes[0].set_title('a.', fontsize=14, loc='left')

# plot foreign population
g = sns.heatmap(afor, annot=True, linewidths=.1, ax=axes[1], cbar=True,
                vmin=0, vmax=1, square=True, cmap='RdPu', fmt='.3f',
                xticklabels=ticks, yticklabels=ticks, cbar_kws={'shrink':0.95})
g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
axes[1].set_title('b.', fontsize=14, loc='left')

# plot somali speakers
g = sns.heatmap(asom, annot=True, linewidths=.1, ax=axes[2], cbar=True,
                vmin=0, vmax=1, square=True, cmap='Oranges', fmt='.3f',
                xticklabels=ticks, yticklabels=ticks, cbar_kws={'shrink':0.95})
g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
axes[2].set_title('c.', fontsize=14, loc='left')

# plot estonian speakers
g = sns.heatmap(aest, annot=True, linewidths=.1, ax=axes[3], cbar=True,
                vmin=0, vmax=1, square=True, cmap='YlGnBu', fmt='.3f',
                xticklabels=ticks, yticklabels=ticks, cbar_kws={'shrink':0.95})
g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
axes[3].set_title('d.', fontsize=14, loc='left')
plt.savefig(r'W:\maphel_langtime\plots\markov_global_comparison_fiswe.pdf', dpi=300,
            bbox_inches='tight')