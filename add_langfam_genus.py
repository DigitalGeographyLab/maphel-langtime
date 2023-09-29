# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 12:50:00 2022

@author: TuoVaisanen-e01
"""
import pandas as pd
import geopandas as gpd

# loop over years
for i in range (1987,2020):
    
    # print messages
    print('[INFO] - Processing year ' + str(i) + '...')
    
    # read hma data in
    hma = gpd.read_file(r'W:\maphel_langtime\geopackage\HMA_langs_div_{}.gpkg'.format(str(i)))
    
    # read national divs in
    fam = pd.read_pickle('W:\\language\\spatial\\pickles\\250m\\harmonized\\' + str(i) + '_langfam_diversity_euref250.pkl')
    gen = pd.read_pickle('W:\\language\\spatial\\pickles\\250m\\harmonized\\' + str(i) + '_langgen_diversity_euref250.pkl')
    
    # join data
    hma = pd.merge(hma, fam, on='euref_250', how='left')
    hma = pd.merge(hma, gen, on='euref_250', how='left')
    
    # save to geopackge
    print('[INFO] - Saving year ' + str(i) + '...')
    hma.to_file(r'W:\maphel_langtime\geopackage\HMA_langs_famgen_div_{}.gpkg'.format(str(i)),
                driver='GPKG')

# print final message
print('[INFO] - ... done!')