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

# function to calculate spatial distribution and density through the years
def calc_dens(data, popcol, ycol):
    
    # copy the dataframe
    df = data.copy()
    
    # calculate area
    df['area_km'] = df.area / 1000000
    
    # get sum of population per year
    sums = df.groupby(ycol)[popcol].sum().rename('pop').reset_index()
    
    # get total area of all inhabited cells per year
    cells = df.groupby(ycol)['area_km'].sum().reset_index()
    
    # merge dataframes
    result = pd.merge(sums, cells, on=ycol)
    
    # calculate density
    result['density'] = result['pop'] / result['area_km']
    
    return result

# columns for weighted diversities
w_cols = ['sompop', 'estpop', 'pop_count']
tags = ['somali', 'estonian', 'hma avg']

# result list
results = []

# loop over years
for i, year in enumerate(list(range(1987,2020))):
    print('[INFO] - Processing year ' + str(year))

    # read files in
    df = gpd.read_file(r'W:\maphel_langtime\geopackage\HMA_langs_famgen_div_{}.gpkg'.format(str(year)))
    
    # assign year for appending
    df['year'] = year
    
    # divide into hma, est and som subset dataframes
    estdf = df[df['estpop'] >= 1]
    somdf = df[df['sompop'] >= 1]
    hma = df[df['pop_count'] >= 1]
    fordf = df[df['foreign_pop'] >= 1]
    findf = df[df['finpop'] >= 1]
    swedf = df[df['swepop'] >= 1]
    
    for n, subset in enumerate([estdf, somdf, fordf, findf, swedf, hma]):
        # check which pop column to use
        if n == 0:
            popcol = 'estpop'
            tag = 'Estonian-inhabited'
        elif n == 1:
            popcol = 'sompop'
            tag = 'Somali-inhabited'
        elif n == 2:
            popcol = 'foreign_pop'
            tag = 'Foreigner-inhabited'
        elif n == 3:
            popcol = 'finpop'
            tag = 'Finnish-inhabited'
        elif n == 4:
            popcol = 'swepop'
            tag = 'Swedish-inhabited'
        elif n == 5:
            popcol = 'pop_count'
            tag = 'HMA average'
        
        # loop over columns
        for col in ['shannon', 'fam_shannon', 'unique_langs', 'unique_fam']:
            
            # calculate series
            score = subset[col] * subset[popcol]
            
            # get sum
            score = score.sum()
            
            # get population sum
            totpop = subset[popcol].sum()
            
            # get weighted average of metric for current year
            score = score / totpop
            
            # create result dataframe
            result = pd.DataFrame(data={'year':[year],'metric':[col],'score':[score],'type':[tag]})
            
            # append to result df list
            results.append(result)

# concatenate results to dataframe
data = pd.concat(results, ignore_index=True)

# save to disk
data.to_pickle(r'W:\maphel_langtime\pickles\popweigh_divs.pkl')
data.to_csv(r'W:\maphel_langtime\pickles\popweigh_divs.csv', sep=';', encoding='utf-8')

# read plotdf
plotdf = pd.read_pickle(r'W:\maphel_langtime\pickles\normalized.pkl')

# set seaborn theme
sns.set()

# plot weighted averages of linguistic diversity
fig, axes = plt.subplots(2,2, figsize=(15,11))
axes = axes.flatten()
palette = {'Somali-inhabited':'C1','Estonian-inhabited':'C0','HMA average':'C2',
           'Finnish-inhabited':'C3','Foreigner-inhabited':'C4', 'Swedish-inhabited':'C5'}
palette2 = {'som':'C1','est':'C0','hma':'C2','fin':'C3','foreign_':'C4', 'swe':'C5'}
g = sns.lineplot(x='year', y='score', hue='type', legend=False, palette=palette,
                 data=data[data['metric'] == 'unique_langs'], ax=axes[0])
g.set(ylabel='Unique languages', xlabel='', ylim=(0,30))
g = sns.lineplot(x='year', y='score', hue='type', legend=False,palette=palette,
                 data=data[data['metric'] == 'unique_fam'], ax=axes[1])
g.set(ylabel='Unique language families', xlabel='', ylim=(0,8))
g = sns.lineplot(x='year', y='score', hue='type', legend=False,palette=palette,
                 data=data[data['metric'] == 'shannon'], ax=axes[2])
g.set(ylabel='Shannon entropy', xlabel='', ylim=(0,1.63))
g = sns.lineplot(x='year', y='score', hue='type', legend=True,palette=palette,
                 data=data[data['metric'] == 'fam_shannon'], ax=axes[3])
g.set(ylabel='Shannon entropy (lang. fam.)', xlabel='', ylim=(0, 1.63))
plt.legend(title='Residential neighbourhood', loc='lower center',
           bbox_to_anchor=[0.20,0.55])
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_popweigh_diversity_87-19.pdf', dpi=300,
            bbox_inches='tight')   
    
 
# set fontsize for title
fs=15
# plot weighted averages of linguistic diversity with normalized diveristies
fig, axes = plt.subplots(2,2, figsize=(15,11))
axes = axes.flatten()
palette = {'Somali-inhabited':'C1','Estonian-inhabited':'C0','HMA average':'C2',
           'Finnish-inhabited':'C3','Foreigner-inhabited':'C4', 'Swedish-inhabited':'C5'}
palette2 = {'som':'C1','est':'C0','hma':'C2','fin':'C3','foreign_':'C4', 'swe':'C5'}
g = sns.lineplot(x='year', y='score', hue='type', legend=False, palette=palette,
                 data=data[data['metric'] == 'unique_langs'], ax=axes[0])
g.set(ylabel='Unique languages', xlabel='', ylim=(0,30))
g.set_title('a.', loc='left', fontsize=fs)
g = sns.lineplot(x='year', y='norm_unique', hue='type', legend=False,
                 palette=palette2, data=plotdf, ci=99, n_boot=3000, ax=axes[1])
g.set(ylabel='Normalized Unique languages', xlabel='', ylim=(-0.3,4))
g.set_title('b.', loc='left', fontsize=fs)
g = sns.lineplot(x='year', y='score', hue='type', legend=False,palette=palette,
                 data=data[data['metric'] == 'shannon'], ax=axes[2])
g.set(ylabel='Shannon entropy', xlabel='', ylim=(0,1.63))
g.set_title('c.', loc='left', fontsize=fs)
g = sns.lineplot(x='year', y='norm_shannon', hue='type', legend=True,
                 palette=palette2, data=plotdf, ci=99, n_boot=3000, ax=axes[3])
g.set(ylabel='Normalized Shannon entropy', xlabel='', ylim=(-0.3,2.5))
g.set_title('d.', loc='left', fontsize=fs)
plt.legend(title='Residential neighbourhood', bbox_to_anchor=[0.80,0.65],
           loc='lower center', labels=['Estonian-inhabited','Somali-inhabited',
                                      'Foreign-inhabited','Finnish-inhabited',
                                      'Swedish-inhabited','HMA average'])
fig.savefig(r'W:\maphel_langtime\plots\popweigh_and_norm_diversity_87-19.pdf', dpi=300,
            bbox_inches='tight')   
    










# plot living space categories
fig, ax = plt.subplots(figsize=(13,10))
palette = {'Somali':'C1','Estonian':'C0','HMA avg.':'C2'}
b = sns.lineplot(data=ls, x='year', y='count', hue='type', style='asva', palette=palette)
b.set(yscale='log', xlabel='', ylabel='Inhabitants')
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_livingspace_lineplot_87-19.pdf', dpi=300,
            bbox_inches='tight')



# plotti missä finsweprop, verrannollinen prop, tu
# earned income is available only through years 1987-2018
# use np.nanmean as estimator to deal with NaNs and confidence intervals
fig, axes = plt.subplots(3, 2, figsize=(15,15))
axes = axes.flatten()
boots = 3000
g = sns.lineplot(x='year', y='finswe_norm', hue='type', data=plotdf, ax=axes[0],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Normalized Finnish/Swedish proportion', xlabel='')
g = sns.lineplot(x='year', y='hi_ed_norm', hue='type', data=plotdf, ax=axes[1],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Normalized high education', xlabel='')
g = sns.lineplot(x='year', y='howner_norm', hue='type', data=plotdf, ax=axes[2],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Normalized owner occupancy rates', xlabel='')
g = sns.lineplot(x='year', y='unemp_norm', hue='type', data=plotdf, ax=axes[3],
             ci=99, n_boot=boots, estimator=np.nanmean,
             legend=False)
g.set(ylabel='Normalized unemployment', xlabel='')
g = sns.lineplot(x='year', y='commute_norm', hue='type', data=plotdf,
                 ax=axes[4], ci=99, n_boot=boots, estimator=np.nanmean,
                 legend=False)
g.set(ylabel='Normalized commute distances', xlabel='')
g = sns.lineplot(x='year', y='income_norm', hue='type', data=plotdf,
                 ax=axes[5], ci=99, n_boot=boots, estimator=np.nanmean,
                 legend=False)
g.set(ylabel='Normalized avg. income', xlabel='')

# define legend and legend location
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'HMA average'],
           bbox_to_anchor=[0.81,1.81])

# save the figure to disk
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_norm_lineplot_87-19.pdf', dpi=300,
            bbox_inches='tight')




# plot density and spatial coverage plot
fig, ax = plt.subplots(1,3, figsize=(19,7))
ax = ax.flatten()
b = sns.lineplot(x='year', y='density', hue='type', estimator=np.nanmean,
                 ci=99, n_boot=boots, data=plotdf[plotdf['type'] != 'hma'], ax=ax[0],
                 legend=False)
b.set(ylabel='Avg. population density per grid', xlabel='')
b = sns.lineplot(x='year', y='area_km', hue='type', data=dens, ax=ax[1],
                 legend=False)
b.set(ylabel='Total area of inhabited grids (km2)', xlabel='')
b = sns.lineplot(x='year', y='finswe_prop', hue='type', estimator=np.nanmean,
                 ci=99, data=plotdf, ax=ax[2], legend=False)
b.set(ylabel='Proportion of Fin/Swe speakers in grid cell', xlabel='')
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'HMA average'],
           bbox_to_anchor=[0.7,0.79])
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_dens_dist_87-19.pdf', dpi=300,
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
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_lineplot_langfam_prop_87-19.pdf', dpi=300,
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
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_lineplot_lf_per_residential_areas_87-19.pdf', dpi=300,
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
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_lineplot_uniques_87-19.pdf', dpi=300,
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
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_lineplot_shannon_87-19.pdf', dpi=300,
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
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_lineplot_diversity_87-19.pdf', dpi=300,
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
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_lineplot_langfam_87.19.pdf',
            dpi=300, bbox_inches='tight')



# plot regression of the linguistic diversity markers
fig, axes = plt.subplots(3, 2, figsize=(14,15))
axes = axes.flatten()

# set up dictionaries for columns and labels
divdict = {0:'unique_langs', 1:'unique_fam', 2:'shannon', 3:'fam_shannon',
           4:'norm_shannon', 5:'norm_fam_shannon'}
ylabels = {0:'Unique languages', 1:'Unique language families', 2:'Shannon entropy (indv. lang.)',
           3:'Shannon entropy (lang. fam.)', 4:'Normalized Shannon entropy (indv. lang.)',
           5:'Normalized Shannon entropy (lang. fam.)'}

# plot regression plots in for loop
for i, ax in enumerate(axes):
    
    # get diversity index
    div = divdict[i]
    
    # loop over dataframes
    for j, data in enumerate([estdf, somdf, hmadf]):
        
        # plot regression
        g = sns.regplot(data=data, x='year', y=div, ci=99, order=3, n_boot=2000,
                        robust=False, scatter=False, ax=ax)
        g.set(ylabel=ylabels[i], xlabel='')
        
# define legend and legend location
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'HMA average'],
           bbox_to_anchor=[0.21,1.85])

# save figure to disk
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_div_regplot_87-19.pdf', dpi=300,
            bbox_inches='tight')





# plot robust regression of the linguistic diversity markers
fig, axes = plt.subplots(3, 2, figsize=(14,15))
axes = axes.flatten()

# set up dictionaries for columns and labels
divdict = {0:'unique_langs', 1:'unique_fam', 2:'shannon', 3:'fam_shannon',
           4:'norm_shannon', 5:'norm_fam_shannon'}
ylabels = {0:'Unique languages', 1:'Unique language families', 2:'Shannon entropy (indv. lang.)',
           3:'Shannon entropy (lang. fam.)', 4:'Normalized Shannon entropy (indv. lang.)',
           5:'Normalized Shannon entropy (lang. fam.)'}

# plot regression plots in for loop
for i, ax in enumerate(axes):
    
    # get diversity index
    div = divdict[i]
    
    # loop over dataframes
    for j, data in enumerate([estdf, somdf, hmadf]):
        
        # plot regression
        g = sns.regplot(data=data, x='year', y=div, ci=99, order=1, n_boot=100,
                        robust=True, scatter=False, ax=ax)
        g.set(ylabel=ylabels[i], xlabel='')
        
# define legend and legend location
plt.legend(title='Residential environment', loc='lower center',
           labels=['Estonian-inhabited', 'Somali-inhabited', 'HMA average'],
           bbox_to_anchor=[0.21,1.85])

# save figure to disk
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_robust_regplot_87-19.pdf', dpi=300,
            bbox_inches='tight')





# plot regression of the socioeconomic markers
fig, axes = plt.subplots(3, 2, figsize=(13,14))
axes = axes.flatten()

# set up dictionaries for columns and labels
divdict = {0:'finswe_prop', 1:'hi_ed_prop', 2:'howner_prop', 3:'unemp_prop',
           4:'avg_earned_income', 5:'income_norm'}
ylabels = {0:'Finnish/Swedish (%)', 1:'Highly educated (%)', 2:'Owner occupancy (%)',
           3:'Unemployment (%)', 4:'Average income (€)', 5:'Normalized avg. income'}

# plot regression plots in for loop
for i, ax in enumerate(axes):
    
    # get diversity index
    div = divdict[i]
    
    # loop over dataframes
    for j, data in enumerate([estdf, somdf, hmadf]):
        
        # ensure everything is float
        data[div] = data[div].astype(float)
        
        # plot regression
        g = sns.regplot(data=data, x='year', y=div, ci=99, order=3,
                        n_boot=2000, robust=False, scatter=False, ax=ax)
        g.set(ylabel=ylabels[i], xlabel='')
        
# define legend and legend location
plt.legend(title='Neighbourhood', loc='lower center',
           labels=['Estonian', 'Somali', 'HMA average'],
           bbox_to_anchor=[0.8,1.9])

# save figure to disk
fig.savefig(r'W:\maphel_langtime\plots\som_est_hma_regplot_socioeco_87-19.pdf', dpi=300,
            bbox_inches='tight')
