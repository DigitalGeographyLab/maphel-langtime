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
import argparse

# set up argument parser
ap = argparse.ArgumentParser()

# Get grid file
ap.add_argument("-g", "--grid", required=True,
                help="Path to folder containing grid history files. For example: /path/to/folder/")

# Get path to input file
ap.add_argument("-sg", "--somagrid", required=True,
                help="Path to folder containing geopackages with diversity metrics calculated. Files should be named 'HMA_langs_famgen_div_[YEAR].gpkg'")

# Get path to input file
ap.add_argument("-eg", "--estogrid", required=True,
                help="Path to folder containing geopackages with diversity metrics calculated. Files should be named 'HMA_langs_famgen_div_[YEAR].gpkg'")

# Get path to input file
ap.add_argument("-fg", "--fingrid", required=True,
                help="Path to folder containing CSVs about commutes.")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output folder. For example: /path/to/folder/. This script assumes you have access to FOLK data within Fiona")

# parse arguments
args = vars(ap.parse_args())

def get_cell_trajectories(dataframe):
    
    # first drop NaNs
    stable = dataframe.dropna()
    
    # get unstable grid cells
    unstable = dataframe[~dataframe['euref_250'].isin(stable['euref_250'].values.tolist())]
    
    # get new and old cells
    new = unstable.dropna(subset=['2019'])
    old = unstable[~unstable['euref_250'].isin(new['euref_250'].values.tolist())]
    
    # return dataframes
    return stable, new, old

# define file
for file in ['shannon', 'unique_langs']:

    # read grid history files in
    df = gpd.read_file(args['grid'] + 'norm_{}_grid_history.gpkg'.format(file))
    sdf = gpd.read_file(args['somagrid'] + 'sompop_grid_history.gpkg')
    edf = gpd.read_file(args['estogrid'] + 'estpop_grid_history.gpkg')
    fdf = gpd.read_file(args['fingrid'] + 'finpop_grid_history.gpkg')

    # set somali annual range to start from 1992 to ensure enough observations
    sdf = sdf.drop(columns=['1987','1988','1989','1990','1991'])

    # drop all rows that have only zero observations throughout
    # to remove cells without any somali or estonian speaking inhabitants
    sdf = sdf[sdf[sdf.columns[1:29]].sum(axis=1) > 0]
    edf = edf[edf[edf.columns[1:34]].sum(axis=1) > 0]
    fdf = fdf[fdf[fdf.columns[1:34]].sum(axis=1) > 0]

    # replace zeros with nans for estonian and somali dataframes
    edf = edf.replace(0, np.nan)
    sdf = sdf.replace(0, np.nan)
    fdf = fdf.replace(0, np.nan)

    # get dataframe triplets
    hma_ey, hma_new, hma_old = get_cell_trajectories(df)
    est_ey, est_new, est_old = get_cell_trajectories(edf)
    som_ey, som_new, som_old = get_cell_trajectories(sdf)
    fin_ey, fin_new, fin_old = get_cell_trajectories(fdf)

    # empty dictionaries for grid ids
    s_eydict = {}
    s_newdict = {}
    s_olddict = {}
    e_eydict = {}
    e_newdict = {}
    e_olddict = {}
    f_eydict = {}
    f_newdict = {}
    f_olddict = {}

    # loop over dataframes
    for i, data in enumerate([som_ey, som_new, som_old, est_ey, est_new, est_old, 
                            fin_ey, fin_new, fin_old]):
        
        # set year range
        if i <= 2:
            yrange = range(1992,2020)
        else:
            yrange = range(1987,2020)
        
        # loop over years
        for year in yrange:
            
            # get string year
            yr = str(year)
            
            # get annual values
            curdf = data[['euref_250', yr]]
            
            # get only inhabited cells
            curdf = curdf[curdf[yr] >= 1]
            
            # get grid cell ids
            gids = curdf['euref_250'].values.tolist()
            
            # append correct dict
            if i == 0:
                s_eydict[year] = gids
            elif i == 1:
                s_newdict[year] = gids
            elif i == 2:
                s_olddict[year] = gids
            elif i == 3:
                e_eydict[year] = gids
            elif i == 4:
                e_newdict[year] = gids
            elif i == 5:
                e_olddict[year] = gids
            elif i == 6:
                f_eydict[year] = gids
            elif i == 7:
                f_newdict[year] = gids
            elif i == 8:
                f_olddict[year] = gids


    # get grid cells with observations for every year
    ey = df.dropna()

    # get grid cells without observations for every year
    wey = df[~df['euref_250'].isin(ey['euref_250'].values.tolist())]

    # get grid cells that are "new" and cells that ceased to exist
    new = wey.dropna(subset=['2019'])
    old = wey[~wey['euref_250'].isin(new['euref_250'].values.tolist())]

    # set list for yearly dfs
    yeardfs = []

    # set type dictionary for yearly difs
    types = {0:'Present every year', 1:'Ceased grid cells', 2:'New grid cells'}

    # loop over dataframes
    for i, data in enumerate([ey, old, new]): 
        
        # define data type
        datatype = types[i]
        
        # loop over years
        for year in range(1987,2020):
            
            # set to string
            yr = str(year)
            
            # get general subset
            subset = data[['euref_250', yr]]
            
            # modify subset for plotting convenience
            subset = subset.rename(columns={yr:file})
            subset['year'] = year
            subset['type'] = datatype
            subset['gridtype'] = 'HMA avg.'
            
            # append to dataframe
            #yeardfs.append(subset)

    # loop over years
    for year in range(1987,2020):
        
        # set to string
        yr = str(year)
        
        # get subset for language specifics
        langsub = df[['euref_250', yr]]
        langsub = langsub.rename(columns={yr:file})
        langsub['year'] = year
        
        # get somali and estonia subsets
        if year < 1992:
            estey = langsub[langsub['euref_250'].isin(e_eydict[year])]
            estey['type'] = 'Present every year'
            estey['gridtype'] = 'Estonian-inhabited'
            estnew = langsub[langsub['euref_250'].isin(e_newdict[year])]
            estnew['type'] = 'New grid cells'
            estnew['gridtype'] = 'Estonian-inhabited'
            estold = langsub[langsub['euref_250'].isin(e_olddict[year])]
            estold['type'] = 'Ceased grid cells'
            estold['gridtype'] = 'Estonian-inhabited'
            # finnish
            finey = langsub[langsub['euref_250'].isin(f_eydict[year])]
            finey['type'] = 'Present every year'
            finey['gridtype'] = 'Finnish-inhabited'
            finnew = langsub[langsub['euref_250'].isin(f_newdict[year])]
            finnew['type'] = 'New grid cells'
            finnew['gridtype'] = 'Finnish-inhabited'
            finold = langsub[langsub['euref_250'].isin(f_olddict[year])]
            finold['type'] = 'Ceased grid cells'
            finold['gridtype'] = 'Finnish-inhabited'
            # append to list
            yeardfs.append(estey)
            yeardfs.append(estnew)
            yeardfs.append(estold)
            yeardfs.append(finey)
            yeardfs.append(finnew)
            yeardfs.append(finold)
        else:
            somey = langsub[langsub['euref_250'].isin(s_eydict[year])]
            somey['type'] = 'Present every year'
            somey['gridtype'] = 'Somali-inhabited'
            somnew = langsub[langsub['euref_250'].isin(s_newdict[year])]
            somnew['type'] = 'New grid cells'
            somnew['gridtype'] = 'Somali-inhabited'
            somold = langsub[langsub['euref_250'].isin(s_olddict[year])]
            somold['type'] = 'Ceased grid cells'
            somold['gridtype'] = 'Somali-inhabited'
            estey = langsub[langsub['euref_250'].isin(e_eydict[year])]
            estey['type'] = 'Present every year'
            estey['gridtype'] = 'Estonian-inhabited'
            estnew = langsub[langsub['euref_250'].isin(e_newdict[year])]
            estnew['type'] = 'New grid cells'
            estnew['gridtype'] = 'Estonian-inhabited'
            estold = langsub[langsub['euref_250'].isin(e_olddict[year])]
            estold['type'] = 'Ceased grid cells'
            estold['gridtype'] = 'Estonian-inhabited'
            # finnish
            finey = langsub[langsub['euref_250'].isin(f_eydict[year])]
            finey['type'] = 'Present every year'
            finey['gridtype'] = 'Finnish-inhabited'
            finnew = langsub[langsub['euref_250'].isin(f_newdict[year])]
            finnew['type'] = 'New grid cells'
            finnew['gridtype'] = 'Finnish-inhabited'
            finold = langsub[langsub['euref_250'].isin(f_olddict[year])]
            finold['type'] = 'Ceased grid cells'
            finold['gridtype'] = 'Finnish-inhabited'
        
            # append to dataframe
            yeardfs.append(somey)
            yeardfs.append(somnew)
            yeardfs.append(somold)
            yeardfs.append(estey)
            yeardfs.append(estnew)
            yeardfs.append(estold)
            yeardfs.append(finey)
            yeardfs.append(finnew)
            yeardfs.append(finold)

    # concatenate dataframes
    result = pd.concat(yeardfs)
    result.to_pickle(args['output'] + 'new_old_neighbourhoods_norm_{}.pkl'.format(file))

    # drop ceased
    result = result[result['type'] != 'Ceased grid cells']

    # plot temporal development of grid cells
    sns.set(font_scale=1.3)
    fig, ax = plt.subplots(figsize=(7,6))
    palette = {'Somali-inhabited':'C1','Estonian-inhabited':'C0','Finnish-inhabited':'C3'}
    g = sns.lineplot(x='year', y='shannon', hue='gridtype', style='type', data=result,
                    ci=99, estimator=np.nanmean, n_boot=3000, palette=palette, ax=ax,
                    legend=False)
    g.set(xlabel='', ylabel='Normalized Shannon entropy')
    #handles = ax.get_legend().legendHandles
    #ax.get_legend().remove()
    #ax.legend(handles, ['Residential environment', 'HMA average',
    #                    'Estonian-inhabited','Somali-inhabited','\nTemporal type',
    #                    'Present every year', 'New grid cells'],
    #          bbox_to_anchor=[0,1.04], loc='upper left')
    plt.savefig(args['output'] + 'grid_trajectory_norm_shannon.pdf', dpi=300,
                bbox_inches='tight')


    # plot temporal development of grid cells unique languages
    fig, ax = plt.subplots(figsize=(7,6))
    palette = {'Somali-inhabited':'C1','Estonian-inhabited':'C0','Finnish-inhabited':'C3'}
    g = sns.lineplot(x='year', y='unique_langs', hue='gridtype', style='type', data=result,
                    ci=99, estimator=np.nanmean, n_boot=3000, palette=palette,
                    legend=False, ax=ax)
    g.set(xlabel='', ylabel='Normalized unique languages')
    plt.savefig(args['output'] + 'grid_trajectory_norm_unique.pdf', dpi=300,
                bbox_inches='tight')