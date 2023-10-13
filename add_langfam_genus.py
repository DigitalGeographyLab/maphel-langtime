# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 12:50:00 2022

@author: TuoVaisanen-e01
"""
import pandas as pd
import geopandas as gpd

import argparse

# set up argument parser
ap = argparse.ArgumentParser()

# Get path to input file
ap.add_argument("-fam", "--family", required=True,
                help="Path to input 250 m grid with calculated diversity based on language families. The filenames need to be '[YEAR]_langfam_diversity_euref250.pkl.")

# Get path to input file
ap.add_argument("-g", "--grid", required=True,
                help="Path to folder containing grid geopackages per year. The filenames need to be 'HMA_langs_div_[YEAR].gpkg")

# Get path to output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output folder. For example: /path/to/folder/. The files are called 'HMA_langs_fam_div_[YEAR].gpkg'")

# parse arguments
args = vars(ap.parse_args())

# get paths
gridpath = args['grid']
fampath = args['family']
outpath = args['output']

# loop over years
for i in range (1987,2020):
    
    # print messages
    print('[INFO] - Processing year ' + str(i) + '...')
    
    # read hma data in
    hma = gpd.read_file(gridpath + 'HMA_langs_div_{}.gpkg'.format(str(i)))
    
    # read national divs in
    fam = pd.read_pickle(fampath + str(i) + '_langfam_diversity_euref250.pkl')
    
    # join data
    hma = pd.merge(hma, fam, on='euref_250', how='left')
    
    # save to geopackge
    print('[INFO] - Saving year ' + str(i) + '...')
    hma.to_file(r'W:\maphel_langtime\geopackage\HMA_langs_famgen_div_{}.gpkg'.format(str(i)),
                driver='GPKG')

# print final message
print('[INFO] - ... done!')