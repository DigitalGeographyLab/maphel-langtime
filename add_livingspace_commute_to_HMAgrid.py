# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 13:22:32 2022

@author: TuoVaisanen-e01
"""
import pandas as pd
import geopandas as gpd
import argparse

# set up argument parser
ap = argparse.ArgumentParser()

# Get grid file
ap.add_argument("-g", "--grid", required=True,
                help="Path to 250 m grid geopackage covering the HMA. For example: /path/to/folder/")

# Get path to input file
ap.add_argument("-lg", "--langgrid", required=True,
                help="Path to folder containing geopackages with diversity metrics calculated. Files should be named 'HMA_langs_famgen_div_[YEAR].gpkg'")

# Get path to input file
ap.add_argument("-hg", "--homegrid", required=True,
                help="Path to folder containing geopackages with diversity metrics calculated. Files should be named 'HMA_langs_famgen_div_[YEAR].gpkg'")

# Get path to input file
ap.add_argument("-c", "--commute", required=True,
                help="Path to folder containing CSVs about commutes.")

# Get path to input file
ap.add_argument("-l", "--livspace", required=True,
                help="Path to folder containing CSVs about living space data.")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output folder. For example: /path/to/folder/. This script assumes you have access to FOLK data within Fiona")

# parse arguments
args = vars(ap.parse_args())

def format_asva(df, dftype, year):
    # calculate types of living arrangements
    df = df['asva'].value_counts().reset_index().rename(columns={'index':'asva','asva':'count'})
    
    # calculate proportions
    df['prop'] = (df['count'] / df['count'].sum()) * 100
    
    # set type and year
    df['type'] = dftype
    df['year'] = year
    
    return df    

# get hma grid
grid = gpd.read_file(args['grid'] + '250m_HMA_accurate.gpkg')

# set grid id to integer
grid['NRO'] = grid['NRO'].astype(int)

# get grid ID's as list
gridlist = list(grid['NRO'].values)

# dataframe list for concatenation
concatlist = []

# loop over years
for i in range(1987,2020):
    print('[INFO] - Processing year ' + str(i))
    
    # get year specific information
    lpath = args['langgrid'] + str(i) + '_mothertongues.csv'
    tpath = args['commute'] + 'FOLK_tkt_data_' + str(i) + '.csv'
    apath = args['livspace'] + 'FOLK_askun_data_' + str(i) + '.csv'
    hpath = args['homegrid'] + 'henkilo_paikkatiedot_' + str(i) +'.csv'
    
    # read data in
    langs = pd.read_csv(lpath, encoding='utf-8', sep=',')
    askun = pd.read_csv(apath, encoding='utf-8', sep=',')
    homes = pd.read_csv(hpath, encoding='utf-8', sep=',')
    
    # drop nan values
    homes = homes.dropna(subset=['euref_250'])
    
    # convert to integer for joining
    homes['euref_250'] = homes['euref_250'].astype(int)
    
    # check if year is not 2019
    if i != 2019:
        # read worklife stats in
        tkt = pd.read_csv(tpath, encoding='utf-8', sep=',')
        # merge tkt and askun data to one
        combined = pd.merge(tkt, askun[['shnro', 'asva']], on='shnro')
        
        # merge language information
        combined = pd.merge(combined, langs[['shnro','kieli']], on='shnro')
    else:
        #merge language information
        combined = pd.merge(askun, langs[['shnro','kieli']], on='shnro')
        
    # merge combined tkt and askun data to homelocation data
    combined = pd.merge(combined, homes[['shnro','euref_250']], on='shnro')
    
    # drop folks who do not live inside the HMA
    HMA = combined[combined['euref_250'].isin(gridlist)]
    
    # get estonians and somalis
    est = HMA[HMA['kieli'] == 'et']
    som = HMA[HMA['kieli'] == 'so']
    
    # get living space per year
    lsest = format_asva(est, 'est', i)
    lssom = format_asva(som, 'som', i)
    lshma = format_asva(HMA, 'hma', i)
    
    # append dataframes
    concatlist.append(lssom)
    concatlist.append(lsest)
    concatlist.append(lshma)
    
    # check if year is not 2019
    if i != 2019:
        
        # get only workforce
        wf = HMA[HMA['ptoim2'] == 11]
        
        # group by grid ids
        print('[INFO] - Grouping by grids for year ' + str(i) + '...')
        grouped = wf.groupby('euref_250')['tyomatka'].mean().reset_index()
    
        
        #save to pickled dataframe
        grouped.to_pickle(args['output'] + 'commute_euclidean_HMA_' + str(i) + '.pkl')
    else:
        pass

# create one df
results = pd.concat(concatlist)
results.to_pickle(args[output] + 'living_space_indv_87-19.pkl')
    
    
    
