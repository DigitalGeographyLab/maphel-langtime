# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:22:03 2022

@author: TuoVaisanen-e01
"""
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def get_cell_trajectories(dataframe):
    
    # first drop NaNs
    stable = dataframe.dropna()
    
    # get unstable grid cells from those not present in stable dataframe
    unstable = dataframe[~dataframe['euref_250'].isin(stable['euref_250'].values.tolist())]
    
    # get new and old cells
    new = unstable.dropna(subset=['2019'])
    old = unstable[~unstable['euref_250'].isin(new['euref_250'].values.tolist())]
    
    # return dataframes
    return stable, new, old

def get_popweigh_values(valdf, popdf, outputtype):
    
    # initialize returnable dataframe
    result = pd.DataFrame()
    
    # loop over years in the data
    for year in range(1987, 2020):
        
        # get year as string
        column = str(year)
        
        # get population sum
        pop = popdf[column].sum()
        
        # get intermediate score
        score = valdf[column] * popdf[column]
        
        # get sum score
        score = score.sum()
        
        # get final score
        score = score / pop
        
        # add result in dataframe
        result.at[year, 'score'] = score
        result.at[year, 'year'] = year
        result.at[year, 'type'] = outputtype
    
    # return the population weighed dataframe
    return result

# initialize a outputlist
outputlist = []

# set labels
labs = ['Unique languages','Shannon entropy']

# loop over files
for x, file in enumerate(['unique_langs','shannon']):

    # read grid history files in
    df = gpd.read_file('W:\\maphel_langtime\\geopackage\\{}_grid_history.gpkg'.format(file))
    sdf = gpd.read_file('W:\\maphel_langtime\\geopackage\\sompop_grid_history.gpkg')
    edf = gpd.read_file('W:\\maphel_langtime\\geopackage\\estpop_grid_history.gpkg')
    hma = gpd.read_file('W:\\maphel_langtime\\geopackage\\pop_count_grid_history.gpkg')
    #finswe = gpd.read_file('W:\\maphel_langtime\\geopackage\\finswe_pop_grid_history.gpkg')
    fin = gpd.read_file('W:\\maphel_langtime\\geopackage\\finpop_grid_history.gpkg')
    #swe = gpd.read_file('W:\\maphel_langtime\\geopackage\\swe_pop_grid_history.gpkg')
    #foreign = gpd.read_file('W:\\maphel_langtime\\geopackage\\foreign_pop_grid_history.gpkg')
    
    # set somali annual range to start from 1992 to ensure enough observations
    sdf = sdf.drop(columns=['1987','1988','1989','1990','1991',])
    #edf = edf.drop(columns=['1987','1988','1989'])
    
    # drop all rows that have only zero observations throughout
    # to remove cells without any somali or estonian speaking inhabitants
    sdf = sdf.drop_duplicates(subset=list(sdf.columns[1:29]), keep=False)
    edf = edf.drop_duplicates(subset=list(edf.columns[1:34]), keep=False)
    fdf = fin.drop_duplicates(subset=list(fin.columns[1:34]), keep=False)
    #swe = swe.drop_duplicates(subset=list(swe.columns[1:31]), keep=False)
    #finswe = finswe.drop_duplicates(subset=list(finswe.columns[1:31]), keep=False)
    #foreign = foreign.drop_duplicates(subset=list(finswe.columns[1:31]), keep=False)
    
    # replace zeros with nans for population count dataframes
    edf = edf.replace(0, np.nan)
    sdf = sdf.replace(0, np.nan)
    fdf = fdf.replace(0, np.nan)
    #swe = swe.replace(0, np.nan)
    #finswe = finswe.replace(0, np.nan)
    #foreign = foreign.replace(0, np.nan)
    
    # get dataframe triplets
    hma_ey, hma_new, hma_old = get_cell_trajectories(hma)
    est_ey, est_new, est_old = get_cell_trajectories(edf)
    som_ey, som_new, som_old = get_cell_trajectories(sdf)
    fin_ey, fin_new, fin_old = get_cell_trajectories(fdf)
    
    # result list
    resultlist = []
    
    # loop over dataframes
    for i, data in enumerate([som_ey, som_new, som_old, est_ey, est_new,
                              est_old, fin_ey, fin_new, fin_old]):
        
        # set year range
        if i <= 2:
            yrange = range(1992,2020)
        elif i > 2:
            yrange = range(1987,2020)
        
        # loop over years
        for year in yrange:
            
            # get current year from pop and value dataframe
            current = data[['euref_250', str(year)]]
            vals = df[['euref_250', str(year)]]
            
            # get grid ids from current year
            gids = current['euref_250'].values.tolist()
            
            # get value grids present in language group dataframe
            vals = vals[vals['euref_250'].isin(gids)]
            
            # get population sum for population weighting
            pop = current[str(year)].sum()
            
            # get intermediate score
            score = vals[str(year)] * current[str(year)]
            
            # get sum score for population weighting
            score = score.sum()
            
            # get final population weighted score
            score = score / pop
            
            # create a dataframe
            result = pd.DataFrame()
            
            # add values to dataframe
            result.at[year, 'year'] = year
            result.at[year, 'score'] = score
            
            # add correct values
            if i == 0:
                result.at[year, 'type'] = 'Somali-inhabited'
                result.at[year, 'celltype'] = 'Present every year'
            elif i == 1:
                result.at[year, 'type'] = 'Somali-inhabited'
                result.at[year, 'celltype'] = 'New'
            elif i == 2:
                result.at[year, 'type'] = 'Somali-inhabited'
                result.at[year, 'celltype'] = 'Ceased'
            elif i == 3:
                result.at[year, 'type'] = 'Estonian-inhabited'
                result.at[year, 'celltype'] = 'Present every year'
            elif i == 4:
                result.at[year, 'type'] = 'Estonian-inhabited'
                result.at[year, 'celltype'] = 'New'
            elif i == 5:
                result.at[year, 'type'] = 'Estonian-inhabited'
                result.at[year, 'celltype'] = 'Ceased'
            elif i == 6:
                result.at[year, 'type'] = 'Finnish-inhabited'
                result.at[year, 'celltype'] = 'Present every year'
            elif i == 7:
                result.at[year, 'type'] = 'Finnish-inhabited'
                result.at[year, 'celltype'] = 'New'
            elif i == 8:
                result.at[year, 'type'] = 'Finnish-inhabited'
                result.at[year, 'celltype'] = 'Ceased'
            
            # append to list
            resultlist.append(result)
    
    # concatenate dataframes
    result = pd.concat(resultlist)
    result.to_pickle(r'W:\maphel_langtime\pickles\new_old_neighbourhoods_popweigh_{}.pkl'.format(file))

    
    # drop rows with ceased
    result = result[result['celltype'] != 'Ceased']
    
    # set name of metric
    result['metric'] = file
    
    # send to output dataframe list
    outputlist.append(result)
    
    # plot temporal development of grid cells
    sns.set(font_scale=1.3)
    fig, ax = plt.subplots(figsize=(7,6))
    palette = {'Somali-inhabited':'C1','Estonian-inhabited':'C0','Finnish-inhabited':'C3'}
    g = sns.lineplot(x='year', y='score', hue='type', style='celltype', data=result,
                     ci=99, estimator=np.nanmean, n_boot=100, palette=palette, ax=ax)
    g.set(xlabel='', ylabel=labs[x])
    if x == 1:
        handles = ax.get_legend().legendHandles
        ax.get_legend().remove()
        ax.legend(handles, ['Residential neighbourhood', 'Somali-inhabited','Estonian-inhabited', 'Finnish-inhabited',
                           '\nTemporal type','Present every year', 'New grid cells'])
    else:
        handles = ax.get_legend().legendHandles
        ax.get_legend().remove()
    plt.savefig(r'W:\maphel_langtime\plots\grid_trajectory_popweigh_{}.pdf'.format(file), dpi=300,
                bbox_inches='tight')

# concatenate the outputs into a dataframe
output = pd.concat(outputlist, ignore_index=True)

# save dataframe
output.to_pickle(r'W:\maphel_langtime\pickles\popweigh_divs_gridtypes.pkl')
output.to_csv(r'W:\maphel_langtime\pickles\popweigh_divs_gridtypes.csv', sep=';',encoding='utf-8')
