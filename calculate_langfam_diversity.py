# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:12:32 2022

@author: TuoVaisanen-e01
"""
import pandas as pd
from collections import Counter
import skbio.diversity.alpha as sk
import numpy as np
import gc
import argparse

# set up argument parser
ap = argparse.ArgumentParser()

# Get path to input file
ap.add_argument("-fam", "--family", required=True,
                help="Path to csv with mappings between language family and languages from Ethnologue")

ap.add_argument("-w", "--wals", required=True,
                help="Path to csv with mappings between language family and languages from WALS")

# Get path to input file
ap.add_argument("-lg", "--langgrid", required=True,
                help="Path to folder containing CSVs with first language information and individual unique ids")

# Get path to input file
ap.add_argument("-hg", "--homegrid", required=True,
                help="Path to folder containing CSVs with home locations and individual unique ids")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output folder. For example: /path/to/folder/. The files will be named: '[YEAR]_langfam_diversity_euref250.pkl'")

# parse arguments
args = vars(ap.parse_args())

# function to count languages in grid and return list of counts
def langcount(langlist):
    count = Counter(langlist)
    return list(count.values())

# read file with mappings between languages and language families
df = pd.read_csv(args['family'])
wals = pd.read_csv(args['wals'])

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

# loop over years 1987-2019 (range ends one before, so 2020 is used)
for i in range(1987,2020):
    
    # get year
    yr = i
    
    # create year-specific filepaths
    lpath = args['langgrid'] + str(yr) + '_mothertongues.csv'
    gpath = args['homegrid'] + 'henkilo_paikkatiedot_' + str(yr) +'.csv'
    outpath = args['outpath'] + str(yr) + '_langfam_diversity_euref250.pkl'
    
    # read annual datasets in
    print('[INFO] - Processing year ' + str(yr) + '...')
    langs = pd.read_csv(lpath, sep=',', encoding='utf-8')
    grids = pd.read_csv(gpath, sep=',', encoding='utf-8')
    
    # merge by shnro
    data = grids.merge(langs[['shnro','kieli']], on='shnro')
    
    # drop langs and grids to release memory
    del langs
    del grids
    
    # count number of people without spatial information
    lost = data['euref_250'].isna().sum()
    print('[INFO] - People lost due to missing spatial data: ' + str(lost))
    
    # drop euref 250 column and NaN rows
    data = data.drop(columns=['euref_1000']).dropna()
    
    # convert grid id to integers
    data['euref_250'] = data['euref_250'].astype(int)
    
    # join language family and genus data
    data = pd.merge(data, df[['alpha2', 'Family', 'Genus']], left_on='kieli', right_on='alpha2', how='left')
    
    # get annual count dataframes
    fam = data['Family'].value_counts().rename(str(yr)).to_frame()
    famdf = pd.concat([famdf, fam], axis=1)
    gen = data['Genus'].value_counts().rename(str(yr)).to_frame()
    gendf = pd.concat([gendf, gen], axis=1)
    
    # group by euref 250 id
    print('[INFO] - Grouping by grids for year ' + str(yr) + '...')
    data = data.groupby('euref_250')['Family'].apply(list).reset_index()
    
    # insert diversity calculations here
    print('[INFO] - Calculating diversities for ' + str(yr) + '...')
    data['pop'] = data['Family'].apply(lambda x: len(x))
    data['unique_fam'] = data['Family'].apply(lambda x: sk.observed_otus(langcount(x)))
    data['fam_berger'] = data['Family'].apply(lambda x: sk.berger_parker_d(langcount(x)))
    data['fam_menhinick'] = data['Family'].apply(lambda x: sk.menhinick(langcount(x)))
    data['fam_shannon'] = data['Family'].apply(lambda x: sk.shannon(langcount(x), base=np.e))
    data['fam_simpson'] = data['Family'].apply(lambda x: sk.simpson(langcount(x)))
    
    # count occurences of

    # drop Genus column
    data = data.drop(columns=['Family'])
    
    # save to pickle
    data.to_pickle(outpath)
    print('[INFO] - Year ' + str(yr) + ' saved to pickle.')
    
    # drop data to release memory and do garbage collect
    del data
    gc.collect()


# plot counts for language families and languages geni across the years
# pivot count datafrmaes
famdf = famdf.stack().reset_index().rename(columns={'level_0':'family','level_1':'year',0:'count'})
gendf = gendf.stack().reset_index().rename(columns={'level_0':'Family','level_1':'year',0:'count'})

# plot
from matplotlib import pyplot as plt
import seaborn as sns

plt.figure(figsize=(15,8))
sns.set(style='whitegrid')
sns.lineplot(data=famdf, x='year', y='count', hue='family', palette='colorblind').set(yscale='log')
plt.xticks(rotation=45)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
plt.savefig(args['outpath'] +'langfam_counts_87-19.pdf')


plt.figure(figsize=(15,8))
sns.set(style='whitegrid')
sns.lineplot(data=gendf, x='year', y='count', hue='genus', palette='colorblind').set(yscale='log')
plt.xticks(rotation=45)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
plt.savefig(args['outpath'] + 'langgen_counts_87-19.pdf')