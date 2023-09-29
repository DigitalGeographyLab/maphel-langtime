# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 12:49:05 2022

@author: TuoVaisanen-e01
"""
import pandas as pd
import gc

# file path list
tkt_paths = ["D:\\ready-made\\FOLK_tkt_8800a\\folk_19872000_tua_tkt21tot_1.dta",
             "D:\\ready-made\\FOLK_tkt_0110a\\folk_20012010_tua_tkt21tot_1.dta",
             "D:\\ready-made\\FOLK_tkt_11a\\folk_20112020_tua_tkt21tot_1.dta",]

ask_paths = ["D:\\ready-made\\FOLK_askun_8800a\\folk_19872000_tua_askun21tot_1.csv",
             "D:\\ready-made\\FOLK_askun_0110a\\folk_20012010_tua_askun21tot_1.csv",
             "D:\\ready-made\\FOLK_askun_11a\\folk_20112019_tua_askun21tot_1.csv"]

# column list to read, tyomatka only available from 2005
tktcols = ['vuosi', 'shnro', 'tyomatka', 'ptoim2']
askcols = ['vuosi', 'shnro', 'asva']

# loop over paths
for path in tkt_paths:
    
    # list for dataframes
    df_list = []
    
    # loop over data and separate annual datasets
    for chunk in pd.read_stata(path, chunksize=100000, columns=tktcols):
        
        # check if year changes in chunk
        if len(list(chunk['vuosi'].unique())) == 2:
            
            # get years present in dataframe
            prevyear = chunk['vuosi'].min()
            newyear = chunk['vuosi'].max()
            
            # print indication message
            print('[INFO] - Reached end of FOLK data for ' + str(prevyear))
            
            # split dataframe in two based on year
            prevdf = chunk[chunk['vuosi'] == prevyear]
            newdf = chunk[chunk['vuosi'] == newyear]
            
            # append previous year to df list
            df_list.append(prevdf)
            
            # concatenate into single dataframe
            data = pd.concat(df_list, ignore_index=True)
            
            # print indication message
            print('[INFO] - Saving FOLK tkt data for ' + str(prevyear) + '...')
            
            # save the full annual data
            data.to_csv('W:\\FOLK\csv\\FOLK_tkt_data_' + str(prevyear) + '.csv',
                        sep=',', encoding='utf-8')
            
            # print indication message
            print('[INFO] - Started processing FOLK tkt data for ' + str(newyear))
            
            # empty dataframe list for next year data
            df_list = []
            
            # add next year data in
            df_list.append(newdf)
            
            # release memory
            del prevdf
            del newdf
            gc.collect()
            
        
        # check if year changes
        elif len(list(chunk['vuosi'].unique())) == 1:
            
            # get current year
            curyear = chunk['vuosi'].max()
            
            
            # append to df list
            df_list.append(chunk)
        
    # concatenate into single dataframe
    data = pd.concat(df_list, ignore_index=True)
    
    # print indication message
    print('[INFO] - Saving FOLK tkt data for ' + str(curyear) + 'e...')
    
    # save
    data.to_csv('W:\\FOLK\\csv\\FOLK_tkt_data_' + str(curyear) + '.csv',
                sep=',', encoding='utf-8')
    
    # release memory
    del data
    gc.collect()

# loop over paths
for path in ask_paths:
    
    # list for dataframes
    df_list = []
    
    # loop over data and separate annual datasets
    for chunk in pd.read_csv(path, sep=',', encoding='utf-8', chunksize=100000,
                             usecols=askcols):
        
        # check if year changes in chunk
        if len(list(chunk['vuosi'].unique())) == 2:
            
            # get years present in dataframe
            prevyear = chunk['vuosi'].min()
            newyear = chunk['vuosi'].max()
            
            # print indication message
            print('[INFO] - Reached end of FOLK askun data for ' + str(prevyear))
            
            # split dataframe in two based on year
            prevdf = chunk[chunk['vuosi'] == prevyear]
            newdf = chunk[chunk['vuosi'] == newyear]
            
            # append previous year to df list
            df_list.append(prevdf)
            
            # concatenate into single dataframe
            data = pd.concat(df_list, ignore_index=True)
            
            # print indication message
            print('[INFO] - Saving FOLK data for ' + str(prevyear) + '...')
            
            # save the full annual data
            data.to_csv('W:\\FOLK\csv\\FOLK_askun_data_' + str(prevyear) + '.csv',
                        sep=',', encoding='utf-8')
            
            # print indication message
            print('[INFO] - Started processing FOLK askun data for ' + str(newyear))
            
            # empty dataframe list for next year data
            df_list = []
            
            # add next year data in
            df_list.append(newdf)
            
            # release memory
            del prevdf
            del newdf
            gc.collect()
            
        
        # check if year changes
        elif len(list(chunk['vuosi'].unique())) == 1:
            
            # get current year
            curyear = chunk['vuosi'].max()
            
            # append to df list
            df_list.append(chunk)
        
    # concatenate into single dataframe
    data = pd.concat(df_list, ignore_index=True)
    
    # print indication message
    print('[INFO] - Saving FOLK data for ' + str(curyear) + 'e...')
    
    # save
    data.to_csv('W:\\FOLK\\csv\\FOLK_askun_data_' + str(curyear) + '.csv',
                sep=',', encoding='utf-8')
    
    # release memory
    del data
    gc.collect()

# print message
print('[INFO] - ... done!')