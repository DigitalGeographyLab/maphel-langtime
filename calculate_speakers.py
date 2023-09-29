# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:58:37 2022

@author: TuoVaisanen-e01
"""
import geopandas as gpd
import pandas as pd
import argparse

# set up argument parser
ap = argparse.ArgumentParser()

# Get path to input file
ap.add_argument("-gp", "--geopackage", required=True,
                help="Path to input 250 m grid (type geopackage).")

# Get path to input file
ap.add_argument("-m", "--mothertongues", required=True,
                help="Path to folder, containing information on first languages in CSVs.")

# Get path to input file
ap.add_argument("-g", "--grid", required=True,
                help="Path to folder containing grid ID information in CSVs.")

# Get path to input file
ap.add_argument("-d", "--diversity", required=True,
                help="Path to folder containing diversity metrics per grid in geopackages.")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output folder. For example: /path/to/folder/")


# parse arguments
args = vars(ap.parse_args())

# get hma grid
hma = gpd.read_file(args['geopackage'])

# get hma grid identifiers
hmagid = list(hma['NRO'].astype(int).values)

# loop over wanted years
for i in range(1987,2020):
    print('[INFO] - Processing year {}...'.format(str(i)))
    
    # create year-specific filepaths
    lpath = args['mothertongues'] + '{}_mothertongues.csv'.format(str(i))
    gpath = args['grid'] + 'henkilo_paikkatiedot_{}.csv'.format(str(i))
    dp = args['diversity'] + 'HMA_folk_processed_language_{}.gpkg'.format(str(i))
    
    # read annual datasets in
    langs = pd.read_csv(lpath, sep=',', encoding='utf-8')
    grids = pd.read_csv(gpath, sep=',', encoding='utf-8')
    
    # simplify grid to hma
    grids = grids.dropna(subset=['euref_250']).drop(columns=['euref_1000'])
    
    # grid id to integers
    grids['euref_250'] = grids['euref_250'].astype(int)
    
    # drop grids outside hma
    grids = grids[grids['euref_250'].isin(hmagid)]
    
    # combine language with home location
    langs = pd.merge(langs[['shnro','kieli']], grids, on='shnro', how='left')
    
    # delet grids from memory
    del grids
    
    # drop non-hma residents
    langs = langs.dropna(subset=['euref_250'])
    
    # grid id to integer
    langs['euref_250'] = langs['euref_250'].astype(int)
    
    # group by grid id
    grouped = langs.groupby(by=['euref_250'])['kieli'].apply(list).reset_index()
    
    # delete langs from memory
    del langs
    
    # calcualate population
    grouped['population'] = grouped['kieli'].apply(len)
    
    # calculate the number of specific speakers
    grouped['finpop'] = grouped['kieli'].apply(lambda x: x.count('fi'))
    grouped['swepop'] = grouped['kieli'].apply(lambda x: x.count('sv'))
    grouped['estpop'] = grouped['kieli'].apply(lambda x: x.count('et'))
    grouped['sompop'] = grouped['kieli'].apply(lambda x: x.count('so'))
    grouped['finswe_pop'] = grouped['finpop'] + grouped['swepop']
    grouped['foreign_pop'] = grouped['population'] - grouped['finswe_pop']
    
    # calculate the proportion of official language speakers
    grouped['finswe_prop'] = round((grouped['finswe_pop'] / grouped['population']) * 100, 3)
    grouped['swe_prop'] = round((grouped['swepop'] / grouped['population']) * 100, 3)
    grouped['fin_prop'] = round((grouped['finpop'] / grouped['population']) * 100, 3) 
    
    # calculate the proportion of the two minority language speakers
    grouped['est_prop'] = round((grouped['estpop'] / grouped['population']) * 100, 3)
    grouped['som_prop'] = round((grouped['sompop'] / grouped['population']) * 100, 3)
    grouped['foreign_prop'] = round((grouped['foreign_pop'] / grouped['population']) * 100, 3)
    
    # calculate concentration of two minority language speakers
    grouped['est_con'] = round((grouped['estpop'] / grouped['estpop'].sum()) * 100, 3)
    grouped['som_con'] = round((grouped['sompop'] / grouped['sompop'].sum()) * 100, 3)
    
    # read diversity grid in
    divgrid = gpd.read_file(dp)
    
    # merge grouped and diversity data
    divgrid = pd.merge(divgrid, grouped[['euref_250', 'estpop', 'sompop',
                                         'finpop', 'swepop', 'finswe_pop',
                                         'foreign_pop', 'finswe_prop', 'est_prop',
                                         'som_prop', 'foreign_prop', 'est_con',
                                         'som_con']],
                       on='euref_250',
                       how='left')
    
    # output path
    o_path = args['output'] + 'HMA_langs_div_{}.gpkg'.format(str(i))
    
    # save geopackage
    print('[INFO] - Saving year {}...'.format(str(i)))
    divgrid.to_file(o_path, driver='GPKG')


print('[INFO] - ... done!')

