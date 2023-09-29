# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:22:14 2022

@author: TuoVaisanen-e01
"""
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns

# set seaborn style
sns.set()

# read file with mappings between languages and language families
mapping = pd.read_csv(r'W:\langfamily_mappings.csv')
wals = pd.read_csv(r'W:\wals_languages.csv')

# simplify data structure by dropping NaN's
mapping = mapping.dropna(subset=['alpha2'])
wals = wals[['ISO639P3code', 'Family', 'Genus']]

# merge
mapping = pd.merge(mapping, wals, left_on='alpha3-b', right_on='ISO639P3code')

# drop duplicates
mapping = mapping.drop_duplicates(subset=['alpha2'])

# simplify
mapping =mapping[['alpha3-b', 'alpha2', 'lang_name', 'ISO639P3code', 'Family', 'Genus']]

# read hma grid data in
hma = gpd.read_file('W:\\grid\\250m_HMA_accurate.gpkg')

# convert grid id to integer
hma['NRO'] = hma['NRO'].astype(int)

# empty lists for results
listlangs = []
listfams = []
listgenus = []

# loop over the years
for i in range(1987, 2020):
    
    # print message
    print('[INFO] - Processing year ' + str(i))
    
    # read data for year
    langs = pd.read_csv('W:\\language\\harmonized\\{}_mothertongues.csv'.format(i),
                        sep=',', encoding='utf-8')
    homes = pd.read_csv('D:\\e01\\custom-made\\henkilo_paikkatiedot_{}.csv'.format(i),
                        sep=',', encoding='utf-8')
    
    homes = homes.dropna(subset=['euref_250'])
    homes['euref_250'] = homes['euref_250'].astype(int)
    langhome = pd.merge(hma[['NRO', 'KUNTA']], homes, how='left', left_on='NRO', right_on='euref_250')
    
    # merge by shnro
    data = langhome.merge(langs[['shnro','kieli']], on='shnro')
    
    # join language family and genus data
    data = pd.merge(data, mapping[['alpha2', 'Family', 'Genus']], left_on='kieli', right_on='alpha2', how='left')
    
    # get counts for languages, language families and genus
    lcounts = data['kieli'].value_counts().rename('count').reset_index().rename(columns={'index':'lang'})
    lfcounts = data['Family'].value_counts().rename('count').reset_index().rename(columns={'index':'lang_fam'})
    gcounts = data['Genus'].value_counts().rename('count').reset_index().rename(columns={'index':'lang_gen'})
    
    # add year information
    lfcounts['year'] = i
    lcounts['year'] = i
    gcounts['year'] = i
    
    # add dataframe to list
    listlangs.append(lcounts)
    listfams.append(lfcounts)
    listgenus.append(gcounts)
    
# concatenate lists to dataframes
res_langs = pd.concat(listlangs)
res_fams = pd.concat(listfams)
res_genus = pd.concat(listgenus)

# group languages by year and get top 10 for every year
toplangs =[]
for year in listlangs:
    top10 = year.head(10)
    toplangs.append(top10)

topgen =[]
for year in listgenus:
    top10 = year.head(10)
    topgen.append(top10)

# concat top 10 languages per year
t10langs = pd.concat(toplangs)
t10gen = pd.concat(topgen)

# plot language development in the HMA
fig, ax = plt.subplots(figsize=(13,10))
g = sns.lineplot(x='year', y='count', hue='lang_fam', data=res_fams)
ax.set(yscale='log')

fig, ax = plt.subplots(figsize=(13,10))
g = sns.lineplot(x='year', y='count', hue='lang', data=t10langs)
ax.set(yscale='log')

fig, ax = plt.subplots(figsize=(13,10))
g = sns.lineplot(x='year', y='count', hue='lang_gen', data=t10gen)
ax.set(yscale='log')