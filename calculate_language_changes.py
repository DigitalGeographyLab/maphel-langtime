# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 15:14:17 2023

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


# set seaborn style
sns.set()

# read hma grid data in
hma = gpd.read_file('W:\\grid\\250m_HMA_accurate.gpkg')

# convert grid id to integer
hma['NRO'] = hma['NRO'].astype(int)

# get hma grid id list
hma_gids = list(hma['NRO'].values)

# list holding dataframes with language information
langlist = []

# loop over the years
for i in range(1987, 2020):
    
    # print message
    print('[INFO] - Processing year ' + str(i))
    
    # read data for year
    langs = pd.read_csv('W:\\language\\harmonized\\{}_mothertongues.csv'.format(i),
                        sep=',', encoding='utf-8')
    homes = pd.read_csv('D:\\e01\\custom-made\\henkilo_paikkatiedot_{}.csv'.format(i),
                        sep=',', encoding='utf-8')
    
    # drop individuals without home grid id
    homes = homes.dropna(subset=['euref_250'])
    
    # convert grid id to integer
    homes['euref_250'] = homes['euref_250'].astype(int)
    
    # get individuals who only live in the HMA
    homes = homes[homes['euref_250'].isin(hma_gids)]
    
    # get languages of individuals who live in the HMA
    langs = pd.merge(homes[['shnro', 'euref_250']], langs[['shnro', 'kieli']],
                     how='left', on='shnro')
    
    # add year
    langs['year'] = i
    
    # save to disk
    langs.to_pickle('W:\\language\\HMA\\HMA_individuals_lang_{}_250m.pkl'.format(str(i)))
    
    # get column new column name
    colname = 'lang' + str(i)
    
    # rename column
    langs = langs.rename(columns={'kieli':colname})
    
    # reindex column to shnro
    langs = langs.set_index('shnro')
    
    # append to list
    langlist.append(langs[colname])
    
# print message
print('[INFO] - Lists appended')

# join dataframes
langhistory = pd.concat(langlist, axis=1)
print('[INFO] - Lists concatenated into a dataframe')

# drop rows with only one obseration
langhistory = langhistory.dropna(thresh=1)

# delete langlist to release memory
del langlist
del hma
del langs
del homes
del hma_gids
gc.collect()

# empty dataframe to hold results
result = pd.DataFrame(columns=['shnro','lang_hist','change_n','change'])

# print message
print('[INFO] - Calculating changes in registered languages....')

# detect changes
for i, row in langhistory.iterrows():
    
    # get language history of individual
    individual = row.dropna().values
    
    # get ordered set of language history (shows language change direction)
    individual = list(dict.fromkeys(individual).keys())
    
    # save only those who have changed
    if len(individual) > 1:
    
        # create one row dataframe
        predf = pd.DataFrame({'shnro':[i],
                              'lang_hist':[individual],
                              'change_n':[len(individual)],
                              'change': [True if len(individual) > 1 else False]})
        
        # add values to dataframe
        result = result.append(predf, ignore_index=True)
    
    else:
        pass

# save resulting dataframe
result.to_pickle('W:\\language\\HMA\\HMA_individuals_lang_changes_1987-2019.pkl')

print('[INFO] - Results saved to pickle!')

# print results
print('\n[INFO] - Top 30 share of language changes from all inhabitants: \n')

print((result['lang_hist'].value_counts()[:30] / len(result)) * 100)

print('\n[INFO] - Top 30 share of final language of all changes: \n')
print((result['lang_hist'].apply(lambda x: x[-1]).value_counts()[:30] / len(result)) * 100)
