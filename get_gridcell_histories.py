# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 12:48:39 2022

@author: TuoVaisanen-e01
"""
import pandas as pd
import geopandas as gpd

# open empty HMA 250 meter grid
grid = gpd.read_file('W:\\grid\\250m_HMA_accurate.gpkg')

# set grid id as integer
grid['euref_250'] = grid['NRO'].astype(int)

# reduce clutter
grid = grid[['euref_250','geometry']]

# loop over relevant column names
for col in ['sompop','estpop','shannon','unique_langs', 'pop_count',
            'finpop', 'swepop','foreign_pop','finswe_pop']:
    
    # set up empty list for dataframes
    dflist = []
    
    # loop over years
    for i in range(1987,2020):
        print('[INFO] - Processing {} data for year {}...'.format(col, str(i)))
    
        # read data
        df = gpd.read_file(r'W:\maphel_langtime\geopackage\HMA_langs_famgen_div_{}.gpkg'.format(i))
        
        # drop uninhabited rows
        df = df[df['pop_count'] >= 1]
        
        # reduce columns to just grid id and the current column from list
        df = df[['euref_250', col]]
        
        # rename current column as the year
        df = df.rename(columns={col:str(i)})
        
        # convert grid id to integer for simpler joining
        df['euref_250'] = df['euref_250'].astype(int)
        
        # append to dataframe list
        dflist.append(df)
    
    # set index for all dataframes as the grid ids
    dflist = [d.set_index('euref_250') for d in dflist]
    
    # concatenate list of dataframes to single dataframe with all the years
    data = pd.concat(dflist, axis=1)
    
    # add geometry information based on the grid ids
    data = pd.merge(data, grid, on='euref_250', how='left')
    
    # convert to geodataframe
    data = gpd.GeoDataFrame(data)
    
    # print message about saving the data to disk
    print('[INFO] - Saving {} column data from 1987 to 2019...'.format(col))
    
    # save to gpkg
    data.to_file('W:\\maphel_langtime\\geopackage\\{}_grid_history.gpkg'.format(col),
                 driver='GPKG')

# Normalized values
# loop over relevant column names
for col in ['sompop','estpop','shannon', 'unique_langs', 'pop_count', 'finpop',
            'swepop','foreign_pop','finswe_pop']:
    
    # set up empty list for dataframes
    dflist = []
    
    # loop over years
    for i in range(1987,2020):
        print('[INFO] - Processing {} data for year {}...'.format(col, str(i)))
    
        # read data
        df = gpd.read_file(r'W:\maphel_langtime\geopackage\HMA_langs_famgen_div_{}.gpkg'.format(i))
        
        # drop uninhabited rows
        df = df[df['pop_count'] >= 1]
        
        # reduce columns to just grid id and the current column from list
        df = df[['euref_250', col]]
        
        # normalize
        df[col] = (df[col] - df[col].mean()) / df[col].std()
        
        # rename current column as the year
        df = df.rename(columns={col:str(i)})
        
        # convert grid id to integer for simpler joining
        df['euref_250'] = df['euref_250'].astype(int)
        
        # append to dataframe list
        dflist.append(df)
    
    # set index for all dataframes as the grid ids
    dflist = [d.set_index('euref_250') for d in dflist]
    
    # concatenate list of dataframes to single dataframe with all the years
    data = pd.concat(dflist, axis=1)
    
    # add geometry information based on the grid ids
    data = pd.merge(data, grid, on='euref_250', how='left')
    
    # convert to geodataframe
    data = gpd.GeoDataFrame(data)
    
    # print message about saving the data to disk
    print('[INFO] - Saving {} column data from 1987 to 2019...'.format(col))
    
    # save to gpkg
    data.to_file('W:\\maphel_langtime\\geopackage\\norm_{}_grid_history.gpkg'.format(col),
                 driver='GPKG')

# print message
print('[INFO] - ... done!')

# get diversities per population group residential neighbourhoods
shannon = gpd.read_file('W:\\maphel_langtime\\geopackage\\shannon_grid_history.gpkg')
somali = gpd.read_file('W:\\maphel_langtime\\geopackage\\sompop_grid_history.gpkg')
estoni = gpd.read_file('W:\\maphel_langtime\\geopackage\\estpop_grid_history.gpkg')
finswe = gpd.read_file('W:\\maphel_langtime\\geopackage\\finswe_pop_grid_history.gpkg')
forgn = gpd.read_file('W:\\maphel_langtime\\geopackage\\foreign_pop_grid_history.gpkg')

# get cells withn population group members
somali = somali.dropna(subset=list(somali.columns[1:34]))
