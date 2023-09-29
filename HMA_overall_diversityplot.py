# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:12:32 2022

@author: TuoVaisanen-e01
"""

import pandas as pd
import geopandas as gpd
from collections import Counter
import skbio.diversity.alpha as sk
import numpy as np
import gc
from matplotlib import pyplot as plt
import seaborn as sns

# function to count languages in grid and return list of counts
def langcount(langlist):
    count = Counter(langlist)
    return list(count.values())

# read file with mappings between languages and language families
print('[INFO] - Reading language mappings...')
df = pd.read_csv(r'W:\langfamily_mappings.csv')
wals = pd.read_csv(r'W:\wals_languages.csv')

# simplify data structure by dropping NaN's
df = df.dropna(subset=['alpha2'])
wals = wals[['ISO639P3code', 'Family', 'Genus']]

# merge
df = pd.merge(df, wals, left_on='alpha3-b', right_on='ISO639P3code')

# drop duplicates
df = df.drop_duplicates(subset=['alpha2'])

# simplify
df = df[['alpha3-b', 'alpha2', 'lang_name', 'ISO639P3code', 'Family', 'Genus']]

# empty dataframes for annual counts of language families and geni
famdf = pd.DataFrame()
gendf = pd.DataFrame()

# read grid in
print('[INFO] - Reading spatial information...')
grid = gpd.read_file('W:\\grid\\250m_HMA_greater_rough.gpkg')

# reduce grid to only HMA
grid = grid[grid['KUNTA'].isin(['091', '049', '092', '235'])]

# convert grid id to integer
grid['NRO'] = grid['NRO'].astype(int)

#get grid ids to list
gids = grid['NRO'].values.tolist()

# result list
resultlist = []

# loop over years 1987-2019 (range ends one before, so 2020 is used)
for i in range(1987,2020):
    print('[INFO] - Processing year {}...'.format(str(i)))
    # get year
    yr = i
    
    # create year-specific filepaths
    lpath = 'W:\\language\\harmonized\\' + str(yr) + '_mothertongues.csv'
    gpath = 'D:\\e01\\custom-made\\henkilo_paikkatiedot_' + str(yr) +'.csv'
    outpath = 'W:\\language\\spatial\\pickles\\250m\\harmonized\\' + str(yr) + '_langfam_diversity_euref250.pkl'
    
    # read annual datasets in
    langs = pd.read_csv(lpath, sep=',', encoding='utf-8')
    grids = pd.read_csv(gpath, sep=',', encoding='utf-8')
    
    # merge on unique individual ids
    data = grids.merge(langs[['shnro','kieli']], on='shnro')
    
    # drop euref 1000 column and then NaN rows in euref 250
    data = data.drop(columns=['euref_1000']).dropna()
    
    # convert grid id to integers
    data['euref_250'] = data['euref_250'].astype(int)
    
    # get only HMA grids
    data = data[data['euref_250'].isin(gids)]
    
    # join language family and genus data
    data = pd.merge(data, df[['alpha2', 'Family', 'Genus']], left_on='kieli', right_on='alpha2', how='left')
    
    # get annual metrics for whole HMA
    shannon = sk.shannon(langcount(data['kieli'].values.tolist()), base=np.e)
    fshannon = sk.shannon(langcount(data['Family'].values.tolist()), base=np.e)
    unique = sk.observed_otus(langcount(data['kieli'].values.tolist()))
    funique = sk.observed_otus(langcount(data['Family'].values.tolist()))
    menhi = sk.menhinick(langcount(data['kieli'].values.tolist()))
    brillouin = sk.brillouin_d(langcount(data['kieli'].values.tolist()))
    heip = sk.heip_e(langcount(data['kieli'].values.tolist()))
    pielou = sk.pielou_e(langcount(data['kieli'].values.tolist()))
    margalef = sk.margalef(langcount(data['kieli'].values.tolist()))
    mcintosh = sk.mcintosh_d(langcount(data['kieli'].values.tolist()))
    simpson = sk.simpson(langcount(data['kieli'].values.tolist()))
    
    fisecount = data['kieli'].value_counts().fi + data['kieli'].value_counts().sv
    focount = len(data) - fisecount
    
    # create dataframe with values
    result = pd.DataFrame()
    
    # add values
    result.at[yr, 'year'] = yr
    result.at[yr, 'shannon'] = shannon
    result.at[yr, 'fam_shannon'] = fshannon
    result.at[yr, 'unique_langs'] = unique
    result.at[yr, 'unique_fam'] = funique
    result.at[yr, 'menhinick'] = menhi
    result.at[yr, 'brillouin'] = brillouin
    result.at[yr, 'heip'] = heip
    result.at[yr, 'pielou'] = pielou
    result.at[yr, 'margalef'] = margalef
    result.at[yr, 'mcintosh'] = mcintosh
    result.at[yr, 'simpson'] = simpson
    result.at[yr, 'natpop'] = fisecount
    result.at[yr, 'forpop'] = focount
    
    # append 
    resultlist.append(result)

# create result dataframe
results = pd.concat(resultlist, ignore_index=True)

# predicted dataframe
preds = {'year':[2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,
                 2032,2033,2034],
         'natpop':[988662,994050,999435,1003592,1006967,1009633,1012031,1014017,
                   1015272,1016151,1016151,1015882,1012678,1008809,1004475],
         'forpop':[216583,228244,240123,252215,264516,277019,289729,302650,
                   315783,329127,342678,356442,370427,384638,399077]}

# create dataframe
preds = pd.DataFrame().from_dict(preds)

# combine preds
results = pd.concat([results, preds], ignore_index=True)

# save results to disk
results.to_pickle(r'W:\maphel_langtime\pickles\HMA_overall_diversities.pkl')
results.to_csv(r'W:\maphel_langtime\pickles\HMA_overall_diversities.csv',
               sep=';', encoding='utf-8')

# initialize seaborn style
sns.set()

# plot the predicted nations
fig, ax = plt.subplots(figsize=(6,5))
g = sns.lineplot(x='year', y='natpop', color='orchid', data=results, ax=ax)
g = sns.lineplot(x='year', y='forpop', color='teal', data=results, ax=ax)
g.set(ylabel='Population', xlabel='')
plt.axvline(2019, 0.05, 0.96, color='k', linestyle='--')
plt.yticks(rotation=45)
plt.ticklabel_format(style='plain', axis='y', useOffset=False)
plt.legend(labels=['Native', 'Non-native'], loc='lower right')
plt.savefig(r'W:\maphel_langtime\plots\hma_population_preds.pdf', dpi=300, bbox_inches='tight')

# plot triple stacked plot for figure 1'
sns.set(font_scale=1.4)
fig, ax = plt.subplots(3,1, figsize=(5,15))
ax = ax.flatten()
g = sns.lineplot(x='year', y='natpop', color='orchid', data=results, ax=ax[0])
g = sns.lineplot(x='year', y='forpop', color='slateblue', data=results, ax=ax[0])
g.set(xlabel='')
g.set_ylabel('Population', fontsize=16)
ax[0].axvline(2019, 0.05, 0.96, color='k', linestyle='--')
ax[0].ticklabel_format(style='plain', axis='y', useOffset=False)
ax[0].tick_params(labelrotation=25)
ax[0].set_yticklabels([-200000, 0, 200000, 400000, 600000, 800000, 1000000, 1200000], rotation=45)
ax[0].legend(labels=['Native', 'Non-native'], loc='center left')
g = sns.lineplot(x='year', y='unique_langs', color='cadetblue', data=results, ax=ax[1])
g.set(xlabel='', ylim=(0,160))
g.set_ylabel('Unique languages', fontsize=16)
g = sns.lineplot(x='year', y='unique_fam', color='olivedrab', data=results, ax=ax[2])
g.set(xlabel='', ylim=(0,14))
g.set_ylabel('Language families', fontsize=16)
plt.savefig(r'W:\maphel_langtime\plots\figure1_bcd.pdf', dpi=300, bbox_inches='tight')
plt.savefig(r'W:\maphel_langtime\plots\figure1_bcd.png', dpi=400, bbox_inches='tight')


# plot global shannon
fig, ax = plt.subplots(1, 2, figsize=(14,5))
ax = ax.flatten()
g = sns.lineplot(x='year', y='shannon', color='olivedrab', data=results, ax=ax[0])
g.set(ylabel='Shannon entropy', xlabel='', ylim=(0,1.3))
g.set_title('a.', loc='left', fontsize=15)
g = sns.lineplot(x='year', y='unique_langs', color='salmon', data=results, ax=ax[1])
g.set(ylabel='Unique languages', xlabel='', ylim=(0,160))
g.set_title('b.', loc='left', fontsize=15)

plt.savefig(r'W:\maphel_langtime\plots\global_shannon_uniques_hma.pdf', dpi=300, bbox_inches='tight')


# plot the global diversities
fig, axes = plt.subplots(2,2, figsize=(12,10))
axes = axes.flatten()
g = sns.lineplot(x='year', y='shannon', data=results, ax=axes[2], color='olivedrab')
g.set(xlabel='', ylabel='Shannon entropy', ylim=(0,1.3))
g = sns.lineplot(x='year', y='fam_shannon', data=results, ax=axes[3], color='teal')
g.set(xlabel='', ylabel='Shannon entropy', ylim=(0,1.3))
g = sns.lineplot(x='year', y='unique_langs', data=results, ax=axes[0], color='salmon')
g.set(xlabel='', ylabel='Languages', ylim=(0,160))
g = sns.lineplot(x='year', y='unique_fam', data=results, ax=axes[1], color='goldenrod')
g.set(xlabel='', ylabel='Language families', ylim=(0,14))
plt.savefig(r'W:\maphel_langtime\plots\afinla_HMA_87-19.pdf', dpi=300, bbox_inches='tight')
