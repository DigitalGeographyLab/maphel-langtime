# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 13:59:54 2022

@author: TuoVaisanen-e01
"""
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import segregation
import argparse

# set up argument parser
ap = argparse.ArgumentParser()

# Get grid file
ap.add_argument("-d", "--diversity", required=True,
                help="Path to folder containing diversities in annual grids. For example: /path/to/folder/")

# Get path to input file
ap.add_argument("-f", "--families", required=True,
                help="Path to folder containing files of language family counts per grid.")

# Get path to input file
ap.add_argument("-c", "--commute", required=True,
                help="Path to folder with annual data on commutes")

# Get path to input file
ap.add_argument("-ls", "--livspace", required=True,
                help="Path to file with living space data between 1987 and 2019.")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output folder. For example: /path/to/folder/. This script assumes you have access to FOLK data within Fiona")

# parse arguments
args = vars(ap.parse_args())

# function to calculate spatial distribution and density through the years
def calc_dens(data, popcol, ycol):
    
    # copy the dataframe
    df = data.copy()
    
    # calculate area
    df['area_km'] = df.area / 1000000
    
    # get sum of population per year
    sums = df.groupby(ycol)[popcol].sum().rename('pop').reset_index()
    
    # get total area of all inhabited cells per year
    cellarea = df.groupby(ycol)['area_km'].sum().reset_index()
    cells = df.groupby(ycol).apply(len).reset_index().rename(columns={0:'cell_count'})
    
    # merge dataframes
    result = pd.merge(sums, cells, on=ycol)
    result = pd.merge(result, cellarea, on=ycol)
    
    # calculate density
    result['density_km'] = result['pop'] / result['area_km']
    result['density_cell'] = result['pop'] / result['cell_count']
    
    return result

# function to calculate relative concentration
def calc_seg(data, yearcol, foccol, popcol):
    
    # copy the initial dataframe
    df = data.copy()
    
    # cast to integers
    df[foccol] = df[foccol].astype(int)
    df[popcol] = df[popcol].astype(int)
    
    # get years from data
    yrmax = df[yearcol].max() + 1
    yrmin = df[yearcol].min()
    
    # result dataframe list
    reslist = []
    
    # loop over years in dataframe
    for year in range(yrmin, yrmax):
        
        # get subset
        subset = df[df[yearcol] == year]
        
        # check whether pop count has discrepancy
        for i, row in subset.iterrows():
            if row[foccol] > row[popcol]:
                print('Fixing discrepancy in population count in {} for {}'.format(str(year), foccol))
                subset.at[i, popcol] = row[foccol]
        
        # calculate aspatial segregation indices
        iso = segregation.aspatial.Isolation(subset, foccol, popcol)
        expo = segregation.aspatial.Exposure(subset, foccol, popcol)
        dissim = segregation.aspatial.Dissim(subset, foccol, popcol)
        atki = segregation.aspatial.Atkinson(subset, foccol, popcol)
        conprof = segregation.aspatial.ConProf(subset, foccol, popcol)
        gini = segregation.aspatial.GiniSeg(subset, foccol, popcol)
        delta = segregation.spatial.Delta(subset, foccol, popcol)
        rco = segregation.spatial.RelativeConcentration(subset, foccol, popcol)
        aco = segregation.spatial.AbsoluteConcentration(subset, foccol, popcol)
        
        
        # get results
        segres = {'isolation':[iso.statistic],'exposure':[expo.statistic],
                  'dissimilarity':[dissim.statistic],'atkinson':[atki.statistic],
                  'concentrationp':[conprof.statistic],'gini':[gini.statistic],
                  'delta':[delta.statistic],'rco':[rco.statistic],
                  'aco':[aco.statistic],'year':[year]}
        
        # add year info
        res = pd.DataFrame().from_dict(data=segres)
        
        # append to list
        reslist.append(res)
    
    # concatenate to single dataframe
    results = pd.concat(reslist)
    
    # return results
    return results

# empty dataframe for catching the annual diversity values from cells with language groups
somlist = []
estlist = []
hmalist = []
finlist = []
swelist = []
forlist = []

# empty list for annual segregation indices
reslist = []

# loop over years
for i in range(1987,2020):
    print('[INFO] - Processing year ' + str(i))
    
    # read files in
    df = gpd.read_file(args['diversity'] +'HMA_langs_famgen_div_{}.gpkg'.format(str(i)))
    famcount = pd.read_pickle(args['families'] + str(i) + '_langfam_counts_HMA_euref250.pkl')
    
    # check if year is ok for commute data
    if i != 2019:
        
        # read pickled commute data
        coms = pd.read_pickle(args['commute'] + 'commute_euclidean_HMA_' + str(i) + '.pkl')
        df = pd.merge(df, coms[['euref_250', 'tyomatka']], on='euref_250', how='left')
    
    # otherwise pass
    else:
        df['tyomatka'] = np.nan
    
    # assign year for appending
    df['year'] = i
    
    # get sum of all work-aged per grid cell
    df['pop_workage'] = df['pop_count'] - (df['pensioners'] + df['under15y'])
    df['total_homes'] = df['home_rental'] + df['home_owner'] + df['home_other'] + df['home_unknown']
    df['total_hi_ed'] = df['ed_bsc'] + df['ed_msc'] + df['ed_phd']
    
    # calculate proportions of ed level, unemployment etc
    df['hi_ed_prop'] = ((df['ed_bsc'] + df['ed_msc'] + df['ed_phd']) / (df['pop_count'] - df['under15y'])) * 100
    df['unemp_prop'] = (df['unemployed'] / df['pop_workage']) * 100
    df['unemp_comp'] = ((df['unemp_prop'] / df['unemp_prop'].mean()) - 1 ) * 100
    df['howner_prop'] = (df['home_owner'] / df['total_homes']) * 100
    df['hrent_prop'] = (df['home_rental'] / df['total_homes']) * 100
    df['income_comp'] = ((df['avg_earned_income'] / df['avg_earned_income'].mean()) - 1) * 100
    df['fin_prop'] = (df['finpop'] / df['pop_count']) * 100
    
    # calculate normalized values
    df['finswe_norm'] = (df['finswe_prop'] - df['finswe_prop'].mean()) / df['finswe_prop'].std()
    df['finprop_norm'] = (df['fin_prop'] - df['fin_prop'].mean()) / df['fin_prop'].std()
    df['forprop_norm'] = (df['foreign_prop'] - df['foreign_prop'].mean()) / df['foreign_prop'].std()
    df['income_norm'] = (df['avg_earned_income'] - df['avg_earned_income'].mean()) / df['avg_earned_income'].std()
    df['entrepre_income_norm'] = (df['avg_entrepre_income'] - df['avg_entrepre_income'].mean()) / df['avg_entrepre_income'].std()
    df['hi_ed_norm'] = (df['hi_ed_prop'] - df['hi_ed_prop'].mean()) / df['hi_ed_prop'].std()
    df['unemp_norm'] = (df['unemployed'] - df['unemployed'].mean()) / df['unemployed'].std()
    df['howner_norm'] = (df['howner_prop'] - df['howner_prop'].mean()) / df['howner_prop'].std()
    df['hrent_norm'] = (df['hrent_prop'] - df['hrent_prop'].mean()) / df['hrent_prop'].std()
    df['commute_norm'] = (df['tyomatka'] - df['tyomatka'].mean()) / df['tyomatka'].std()
    
    # calculate normalize shannon
    # normalize shannon with division by the logarithm of unique observations
    df['norm_shannon'] = (df['shannon'] - df['shannon'].mean()) / df['shannon'].std()
    df['norm_gen_shannon'] = (df['gen_shannon'] - df['gen_shannon'].mean()) / df['gen_shannon'].std()
    df['norm_fam_shannon'] = (df['fam_shannon'] - df['fam_shannon'].mean()) / df['fam_shannon'].std()
    df['norm_unique'] = (df['unique_langs'] - df['unique_langs'].mean()) / df['unique_langs'].std()
    df['norm_fam_unique'] = (df['unique_fam'] - df['unique_fam'].mean()) / df['unique_fam'].std()
    
    # join language family counts
    df = pd.merge(df, famcount, on='euref_250', how='left')

    # get hma dataframe
    hma = df[df['pop_count'] >= 1]
    
    # assign marker for hma
    hma['type'] = 'hma'
    
    # append hma df to list
    hmalist.append(hma)
    
    # get mean proportion of finnish and swedish speakers
    mprop = hma['finswe_prop'].mean()
    
    # loop over languages
    for lang in ['est','som','fin','swe','foreign_']:
        
        # define column names
        if lang != 'pop_count':
            colname = lang + 'pop'
        else:
            colname = 'pop_count'
        
        # get copy of dataframe
        data = hma.copy()
                
        # fill NaNs with 0
        data[colname] = data[colname].fillna(value=0)
        
        # check whether pop count has discrepancy
        for j, row in data.iterrows():
            if row[colname] > row['pop_count']:
                print('Fixing discrepancy in population count in {} for {}'.format(str(i), colname))
                data.at[j, 'pop_count'] = row[colname]
        
        # calculate language-specific segregation metrics
        iso = segregation.aspatial.Isolation(data, colname, 'pop_count').statistic
        expo = segregation.aspatial.Exposure(data, colname, 'pop_count').statistic
        dissim = segregation.aspatial.Dissim(data, colname, 'pop_count').statistic
        atki = segregation.aspatial.Atkinson(data, colname, 'pop_count').statistic
        conprof = segregation.aspatial.ConProf(data, colname, 'pop_count').statistic
        gini = segregation.aspatial.GiniSeg(data, colname, 'pop_count').statistic
        delta = segregation.spatial.Delta(data, colname, 'pop_count').statistic
        
        # get results
        segres = {'isolation':[iso],'exposure':[expo], 'dissimilarity':[dissim],
                  'atkinson':[atki],'concentrationp':[conprof],'gini':[gini],
                  'delta':[delta],'year':[i], 'type':[lang]}
        
        # add year info
        res = pd.DataFrame().from_dict(data=segres)
        
        # append to list
        reslist.append(res)
        
        # get subset
        subset = df[df[colname] >= 1]
        
        # assign markers for subset
        subset['type'] = lang
        
        # calculate comparative proportion
        subset['comprop'] = ((subset['finswe_prop'] / mprop) - 1) * 100
        
        # check which columns to drop, BUT WHY ARE WE DROPPING THEM?
        if lang == 'est':
            
            # append to dataframe
            estlist.append(subset)
            
        elif lang == 'som':
            
            # append to dataframe
            somlist.append(subset)
        
        elif lang == 'fin':
            
            # append to dataframe
            finlist.append(subset)
        
        elif lang == 'swe':
            
            
            # append to dataframe
            swelist.append(subset)
        
        elif lang == 'foreign_':
            
            
            # append to dataframe
            forlist.append(subset)

# create dataframes from list of dataframes
somdf = pd.concat(somlist)
estdf = pd.concat(estlist)
findf = pd.concat(finlist)
swedf = pd.concat(swelist)
fordf = pd.concat(forlist)
hmadf = pd.concat(hmalist)

# actual segregation indices
segre = pd.concat(reslist, ignore_index=True)

# calculate grid densities by somali and estonian speakers
somdf['density'] = somdf['sompop'] / (somdf.area / 1000000)
findf['density'] = findf['finpop'] / (findf.area / 1000000)
swedf['density'] = swedf['swepop'] / (swedf.area / 1000000)
fordf['density'] = fordf['foreign_pop'] / (fordf.area / 1000000)
estdf['density'] = estdf['estpop'] / (estdf.area / 1000000)
hmadf['density'] = hmadf['pop_count'] / (hmadf.area / 1000000)

# create plotdf
plotlist = [estdf, somdf, fordf, findf, swedf, hmadf]
plotdf = pd.concat(plotlist)
plotdf = plotdf.reset_index(drop=True)

# convert income dtype from 'object' to 'float'
plotdf['income_comp'] = plotdf['income_comp'].astype(float)
plotdf['income_norm'] = plotdf['income_norm'].astype(float)
plotdf['avg_earned_income'] = plotdf['avg_earned_income'].astype(float)

# save plotdf for convenience
plotdf.to_pickle(args['output'] + 'language_groups_ready_to_plot.pkl')

# actual segregation indices
segre = pd.concat(reslist, ignore_index=True)

# save segregation to disk
segre.to_pickle(args['output'] + 'segregation.pkl')

# get delta baselines from 1987
baselines = {'est':segre['delta'][0],'som':segre['delta'][1],
             'fin':segre['delta'][2],'swe':segre['delta'][3],'foreign_':segre['delta'][4]}

# calculate delta concentration change from 1987 baseline
for i, row in segre.iterrows():
    segre.at[i, 'delta_comp'] = ((row['delta'] / baselines[row['type']]) - 1) * 100

# calculate spatial distribution of somali and estonian speakers
somdens = calc_dens(somdf, 'sompop', 'year')
somdens['type'] = 'som'
estdens = calc_dens(estdf, 'estpop', 'year')
estdens['type'] = 'est'
hmadens = calc_dens(hmadf, 'pop_count', 'year')
hmadens['type'] = 'hma'
findens = calc_dens(findf, 'finpop', 'year')
findens['type'] = 'fin'
swedens = calc_dens(swedf, 'swepop', 'year')
swedens['type'] = 'swe'
fordens = calc_dens(fordf, 'foreign_pop', 'year')
fordens['type'] = 'foreign_'

# combine together
dens = pd.concat([estdens, somdens, fordens, findens, swedens, hmadens])
dens = dens.reset_index(drop=True)

# save to disk
dens.to_pickle(args['output'] + 'language_group_densities.pkl')

# read in annual living space category count data
ls = pd.read_pickle(args['livspace'] + 'living_space_indv_87-19.pkl')

# dictionary for classifications
asdict = {1:'Spacious',2:'Normal',3:'Overcrowded',4:'Unknown'}
ldict = {'hma':'HMA avg.', 'som':'Somali', 'est':'Estonian'}

# rename
ls = ls.replace({'asva':asdict, 'type':ldict})

# set seaborn theme
sns.set()

# plot living space categories
fig, ax = plt.subplots(figsize=(13,10))
palette = {'Somali':'C1','Estonian':'C0','HMA avg.':'C2'}
b = sns.lineplot(data=ls, x='year', y='count', hue='type', style='asva', palette=palette)
b.set(yscale='log', xlabel='', ylabel='Inhabitants')
fig.savefig(args['output'] + 'som_est_hma_livingspace_lineplot_87-19.pdf', dpi=300,
            bbox_inches='tight')

# use the shared palette
palette = {'som':'C1','est':'C0','hma':'C2',
           'fin':'C3','foreign_':'C4', 'swe':'C5'}
palette2 = {'Somali-inhabited':'C1','Estonian-inhabited':'C0','HMA average':'C2',
           'Finnish-inhabited':'C3','Foreigner-inhabited':'C4', 'Swedish-inhabited':'C5'}
# plotti missä finsweprop, verrannollinen prop, tu
# earned income is available only through years 1987-2018
# use np.nanmean as estimator to deal with NaNs and confidence intervals
fig, axes = plt.subplots(3, 2, figsize=(15,15))
axes = axes.flatten()
boots = 3000
g = sns.lineplot(x='year', y='finprop_norm', hue='type', data=plotdf, ax=axes[0],
             ci=99, n_boot=boots, estimator=np.nanmean, palette=palette,
             legend=False)
g.set(ylabel='Normalized Finnish-speaking proportion', xlabel='')
g.set_title(label='a.', loc='left')
g = sns.lineplot(x='year', y='hi_ed_norm', hue='type', data=plotdf, ax=axes[1],
             ci=99, n_boot=boots, estimator=np.nanmean,palette=palette,
             legend=False)
g.set(ylabel='Normalized high education', xlabel='')
g.set_title(label='b.', loc='left')
g = sns.lineplot(x='year', y='howner_norm', hue='type', data=plotdf, ax=axes[2],
             ci=99, n_boot=boots, estimator=np.nanmean,palette=palette,
             legend=False)
g.set(ylabel='Normalized owner occupancy rates', xlabel='')
g.set_title(label='c.', loc='left')
g = sns.lineplot(x='year', y='income_norm', hue='type', data=plotdf,
                 ax=axes[3], ci=99, n_boot=boots, estimator=np.nanmean,
                 palette=palette,legend=False)
g.set(ylabel='Normalized avg. income', xlabel='')
g.set_title(label='d.', loc='left')
g = sns.lineplot(x='year', y='commute_norm', hue='type', data=plotdf,
                 ax=axes[4], ci=99, n_boot=boots, estimator=np.nanmean,
                 palette=palette,legend=False)
g.set(ylabel='Normalized commute distances', xlabel='')
g.set_title(label='e.', loc='left')
g = sns.lineplot(x='year', y='unemp_norm', hue='type', data=plotdf, ax=axes[5],
             ci=99, n_boot=boots, estimator=np.nanmean,
             palette=palette,legend=False)
g.set(ylabel='Normalized unemployment', xlabel='')
g.set_title(label='f.', loc='left')


# define legend and legend location
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'Foreign-inhabited',
                   'Finnish-inhabited', 'Swedish-inhabited', 'HMA average'],
           bbox_to_anchor=[0.81,0.56])

# save the figure to disk
fig.savefig(args['output'] + 'som_est_hma_norm_lineplot_87-19.pdf', dpi=300,
            bbox_inches='tight')

# get saveable dataframe
savedf = plotdf[['year', 'type', 'finswe_norm', 'finprop_norm', 'forprop_norm',
                 'income_norm', 'entrepre_income_norm', 'hi_ed_norm', 'unemp_norm',
                 'howner_norm', 'hrent_norm', 'commute_norm', 'shannon',
                 'unique_langs', 'norm_shannon', 'norm_fam_shannon',
                 'norm_unique', 'norm_fam_unique']]
savedf = savedf.round(4)
savedf.to_pickle(args['output'] + 'language_groups_normalized_diversities.pkl')

# plot density and spatial coverage plot
fig, ax = plt.subplots(2,2, figsize=(14,14))
ax = ax.flatten()
b = sns.lineplot(x='year', y='density_cell', hue='type', estimator=np.nanmean,
                 ci=99, n_boot=boots, data=dens, ax=ax[0],
                 legend=False, palette=palette)
b.set(ylabel='Avg. population count per grid cell', xlabel='')
b.set_title(label='a.', loc='left')
b = sns.lineplot(x='year', y='delta', hue='type', estimator=np.nanmean,
                 ci=99, n_boot=boots, data=segre, ax=ax[1], legend=False,
                 palette=palette)
b.set(ylabel='Delta concentration index', xlabel='', ylim=(0.37,1.03))
b.set_title(label='b.', loc='left')
b = sns.lineplot(x='year', y='area_km', hue='type', data=dens, ax=ax[2],
                 legend=False, palette=palette)
b.set(ylabel='Total area of inhabited grids (km\u00b2)', xlabel='')
b.set_title(label='c.', loc='left')
b = sns.lineplot(x='year', y='finswe_prop', hue='type', estimator=np.nanmean,
                 ci=99, data=plotdf, ax=ax[3], legend=False, palette=palette)
b.set(ylabel='Proportion of Fin/Swe speakers in grid cell', xlabel='')
b.set_title(label='d.', loc='left')
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'Foreign-inhabited',
                   'Finnish-inhabited','Swedish-inhabited','HMA average'],
           bbox_to_anchor=[0.6,0.1])
fig.savefig(args['output'] + 'som_est_hma_dens_grid_dist_87-19.pdf', dpi=300,
            bbox_inches='tight')




# plot surrounding language families
fig, axes = plt.subplots(2, 2, figsize=(15,13))
axes = axes.flatten()
boots = 3000
g = sns.lineplot(x='year', y='prop_indoeur', hue='type', data=plotdf, ax=axes[0],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Indo-european languages (%)', xlabel='')
g = sns.lineplot(x='year', y='prop_altaic', hue='type', data=plotdf, ax=axes[1],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Altaic languages (%)', xlabel='')
g = sns.lineplot(x='year', y='prop_austroas', hue='type', data=plotdf, ax=axes[2],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Austro-asiatic languages (%)', xlabel='', ylim=(-0.5,7))
g = sns.lineplot(x='year', y='prop_afroas', hue='type', data=plotdf, ax=axes[3],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Afro-asiatic languages (%)', xlabel='')
# define legend and legend location
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'HMA average'],
           bbox_to_anchor=[0.35,0.75])
fig.savefig(args['output'] + 'som_est_hma_lineplot_langfam_prop_87-19.pdf', dpi=300,
            bbox_inches='tight')




# plot language family proportions across somali and estonian grid cells
fig, axes = plt.subplots(1, 2, figsize=(16,7))
axes = axes.flatten()
boots = 3000
titles = ['Grid cells with Somali speakers', 'Grid cells with Estonian speakers']
for i, dataframe in enumerate([somdf, estdf]):
    for j, column in enumerate(['prop_altaic', 'prop_austroas',
                                'prop_afroas', 'prop_taikadai', 'prop_austrones',
                                'prop_nigercong', 'prop_dravidian']):
        g = sns.lineplot(x='year', y=column, data=dataframe, ax=axes[i], ci=99,
                         n_boot=boots, estimator=np.nanmean, legend=False)
        g.set(ylabel='Proportion of inhabitants (%)', xlabel='', ylim=(-0.3,4.5))
plt.legend(title='Language family', loc='lower center',
           labels=['Altaic','Austro-asiatic','Afro-asiatic',
                   'Tai-Kadai','Austronesian', 'Niger-Congo','Dravidian'],
           bbox_to_anchor=(0.5,0.7))
fig.savefig(args['output'] + 'som_est_hma_lineplot_lf_per_residential_areas_87-19.pdf', dpi=300,
            bbox_inches='tight')




# plot unique langs and families absolute and normalized
fig, axes = plt.subplots(2, 2, figsize=(15,13))
axes = axes.flatten()
boots = 3000
g = sns.lineplot(x='year', y='unique_langs', hue='type', data=plotdf, ax=axes[0],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Unique languages', xlabel='', ylim=(0,25))
g = sns.lineplot(x='year', y='unique_fam', hue='type', data=plotdf, ax=axes[1],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Unique language families', xlabel='', ylim=(0,8))
g = sns.lineplot(x='year', y='norm_unique', hue='type', data=plotdf, ax=axes[2],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Unique languages (normalized)', xlabel='', ylim=(-0.3,8))
g = sns.lineplot(x='year', y='norm_fam_unique', hue='type', data=plotdf, ax=axes[3],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Unique language families (normalized)', xlabel='', ylim=(-0.3,6.5))
# define legend and legend location
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'HMA average'],
           bbox_to_anchor=[0.70,0.75])
fig.savefig(args['output'] + 'som_est_hma_lineplot_uniques_87-19.pdf', dpi=300,
            bbox_inches='tight')



# plot shannon langs and families absolute and normalized
fig, axes = plt.subplots(2, 2, figsize=(15,13))
axes = axes.flatten()
boots = 3000
g = sns.lineplot(x='year', y='shannon', hue='type', data=plotdf, ax=axes[0],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Shannon entropy of languages', xlabel='')
g = sns.lineplot(x='year', y='fam_shannon', hue='type', data=plotdf, ax=axes[1],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Shannon entropy of language families', xlabel='')
g = sns.lineplot(x='year', y='norm_shannon', hue='type', data=plotdf, ax=axes[2],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Shannon entropy of languages (normalized)', xlabel='')
g = sns.lineplot(x='year', y='norm_fam_shannon', hue='type', data=plotdf, ax=axes[3],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Shannon entropy of language families (normalized)', xlabel='')
# define legend and legend location
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'HMA average'],
           bbox_to_anchor=[0.70,0.75])
fig.savefig(args['output'] + 'som_est_hma_lineplot_shannon_87-19.pdf', dpi=300,
            bbox_inches='tight')



# plot linguistic diversity with lineplots
fig, axes = plt.subplots(3, 2, figsize=(15,15))
axes = axes.flatten()
boots = 3000
g = sns.lineplot(x='year', y='unique_langs', hue='type', data=plotdf, ax=axes[0],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Unique languages', xlabel='')
g = sns.lineplot(x='year', y='unique_fam', hue='type', data=plotdf, ax=axes[1],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Unique language families', xlabel='')
g = sns.lineplot(x='year', y='shannon', hue='type', data=plotdf, ax=axes[2],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Shannon entropy of languages', xlabel='')
g = sns.lineplot(x='year', y='fam_shannon', hue='type', data=plotdf, ax=axes[3],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Shannon entropy of language families', xlabel='')
g = sns.lineplot(x='year', y='norm_shannon', hue='type', data=plotdf,
                 ax=axes[4], ci=99, n_boot=boots, estimator=np.nanmean,
                 legend=False)
g.set(ylabel='Norm. Shannon entropy', xlabel='', ylim=(-0.5,2.3))
g = sns.lineplot(x='year', y='norm_fam_shannon', hue='type', data=plotdf,
                 ax=axes[5], ci=99, n_boot=boots, estimator=np.nanmean,
                 legend=False)
g.set(ylabel='Norm. Shannon entropy (lf)', xlabel='', ylim=(-0.5,2))

# define legend and legend location
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'HMA average'],
           bbox_to_anchor=[0.70,1.9])
fig.savefig(args['output'] + 'som_est_hma_lineplot_diversity_87-19.pdf', dpi=300,
            bbox_inches='tight')




# plot language families across residential envrionmentss
fig, axes = plt.subplots(1, 3, figsize=(17,7))
axes = axes.flatten()

# set subplots and enumerate
for i, graph in enumerate([somdf,estdf,hmadf]):
    # loop over count columns
    for count in ['prop_indoeur', 'prop_altaic',
                  'prop_austroas', 'prop_afroas', 'prop_taikadai',
                  'prop_austrones', 'prop_nigercong', 'prop_dravidian']:
        # plot counts on top of each other
        g = sns.lineplot(x='year', y=count, data=graph[graph[count] > 0],
                 ax=axes[i], ci=99, n_boot=boots, estimator=np.nanmean,
                 legend=True)
        
    # check which df for title
    if i == 0:
        title = 'Residential environment of Somali speakers'
    elif i == 1:
        title = 'Residential environment of Estonian speakers'
    elif i == 2:
        title = 'HMA average'
    g.set(xlabel='', ylabel='Proportion of inhabitants (%)', title=title,
          ylim=(-0.2,19.5))
    
plt.legend(title='Language family', loc='lower center',
           labels=['Indo-European', 'Altaic', 'Austro-Asiatic', 'Afro-Asiatic',
                   'Tai-Kadai', 'Austronesian', 'Niger-Congo', 'Dravidian'],
           bbox_to_anchor=[0.5,0.3])
fig.savefig(args['output'] + 'som_est_hma_lineplot_langfam_87.19.pdf',
            dpi=300, bbox_inches='tight')