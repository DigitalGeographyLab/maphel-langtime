# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 13:56:11 2022

@author: TuoVaisanen-e01
"""
import giddy
import geopandas as gpd
import pandas as pd
import mapclassify as mc
from libpysal.weights import KNN
import numpy as np
import matplotlib.pyplot as plt
from esda.moran import Moran
import seaborn as sns

# list of files to loop over for normalized values
flist = ['shannon', 'unique_langs', 'finpop', 'swepop',
         'foreign_pop', 'estpop', 'sompop']
tlist = ['Shannon entropy','Unique languages',
         'Finnish speakers', 'Swedish speakers',
         'Foreign-language speakers','Estonian speakers','Somali speakers']
cpalettes = ['PuBuGn','GnBu','OrRd','YlOrBr','RdPu','YlGnBu','YlGn']

# list of files to loop over for normalized values
flist = ['shannon', 'unique_langs']
tlist = ['Shannon entropy','Unique languages']
cpalettes = ['PuBuGn','GnBu']

# variable for quantiles
k = 5

# loop over files
for j, file in enumerate(flist):
    
    # get color palette
    cpal = cpalettes[j]
    
    # print message
    print('[INFO] - Reading in {}...'.format(file))
    
    # read data in
    df = gpd.read_file('W:\\maphel_langtime\\geopackage\\{}_grid_history.gpkg'.format(file))
    
    # get initial spatial weight matrix for the grid KNN=8
    print('[INFO] - Calculating initial KNN-8 weights matrix for Global Moran'\
          ' for {}...'.format(file))
    iw = KNN.from_dataframe(df, k=8)
    
    # get initial cell values
    ivalues = np.array([df[str(col)].values for col in range(1987,2020)])
    ivals = np.nan_to_num(ivalues)
    
    # calculate global moran i per year
    moran = [Moran(year, iw) for year in ivals]
    
    # get results from global morans i
    res = np.array([(mi.I, mi.EI, mi.seI_norm, mi.sim[974]) for mi in moran])
    
    # set initial years
    yrs = np.arange(1987,2020)
    
    # plot moran's I
    fig, ax = plt.subplots(1,1, figsize=(6,5))
    ax.plot(yrs, res[:,0], label="Moran's I")
    ax.plot(yrs, res[:,1] + 1.96 * res[:,2], label='Upper bound', linestyle='dashed')
    ax.plot(yrs, res[:,1] - 1.96 * res[:,2], label='Lower bound', linestyle='dashed')
    ax.set_title("Global Moran's I: " + tlist[j])
    ax.legend()
    plt.savefig(r'W:\maphel_langtime\plots\global_moran_{}.pdf'.format(file),
                dpi=300, bbox_inches='tight')
    
    # get array of data
    if file == 'sompop':
        # drop all rows without any somalis ever
        df = df.loc[(df[df.columns[6:34]].sum(axis=1) != 0)]
        # get values from those rows
        values = np.array([df[str(col)].values for col in range(1992,2020)])
    elif file == 'estpop':
        # do the same for estonians
        df = df.loc[(df[df.columns[6:34]].sum(axis=1) != 0)]
        values = np.array([df[str(col)].values for col in range(1992,2020)])
    elif file == 'foreign_pop':
        # do the same for estonians
        df = df.loc[(df[df.columns[1:34]].sum(axis=1) != 0)]
        values = np.array([df[str(col)].values for col in range(1987,2020)])
    elif file == 'swepop':
        # do the same for estonians
        df = df.loc[(df[df.columns[1:34]].sum(axis=1) != 0)]
        values = np.array([df[str(col)].values for col in range(1987,2020)])
    else:
        df = df.loc[(df[df.columns[1:34]].sum(axis=1) != 0)]
        values = np.array([df[str(col)].values for col in range(1987,2020)])
    
    # fill nan with zeroes
    vals = np.nan_to_num(values)
    
    # set year to test classification scheme
    #yr = 3
    
    # get mapclassify objects
    #quint = mc.Quantiles(vals[yr], k=k)
    #fj = mc.FisherJenks(vals[yr], k=k)
    #ei = mc.EqualInterval(vals[yr], k=k)
    #nb = mc.NaturalBreaks(vals[yr], k=k)
    #mb = mc.MaximumBreaks(vals[yr], k=k)
    #bp = mc.BoxPlot(vals[yr])
    #msd = mc.StdMean(vals[yr])
    
    # group classifier objects
    #classes = quint, ei, fj, nb, mb, bp, msd
    
    # get adcm for classifiers
    #fits = np.array([c.adcm for c in classes])
    
    # convert to df
    #adcms = pd.DataFrame(fits)
    
    # add names
    #adcms['classifier'] = [c.name for c in classes]
    
    # set colujmn names
    #adcms.columns = ['ADCM', 'Classifier']
    
    # plot barplot
    #fig, ax = plt.subplots(figsize=(10,8))
    #g = sns.barplot(y='Classifier', x='ADCM', data=adcms, palette='Pastel1')
    #g.set(title=tlist[j])
    #plt.savefig(r'W:\maphel_langtime\plots\spatial_markov_adcms'\
    #            '{}_k{}.pdf'.format(file,str(k)),
    #            dpi=300, bbox_inches='tight')
        
    # NB: JenksCaspall is the best classifier, but breaks down in this for some reason
    # Second best is FisherJenks or NaturalBreaks, but Quantiles is the easiest
    # to understand
    # classify values into quintiles per year
    print('[INFO] - Classifying every year for {}...'.format(file))
    #q5 = np.array([mc.Quantiles(year, k=k).yb for year in vals]).transpose()
    #q5 = np.array([mc.Percentiles(year, pct=[20,40,60,80,100]).yb for year in vals]).transpose()
    #q5 = np.array([mc.FisherJenks(year, k=k).yb for year in vals]).transpose()
    q5 = np.array([mc.NaturalBreaks(year, k=k, initial=40).yb for year in vals]).transpose()
    #q5 = np.array([mc.JenksCaspall(year, k=k).yb for year in vals]).transpose()
    # Use Natural Breaks for NORMALIZED VALUES! QUINTILES DO NOT WORK!
    
    # list to append dataframe
    classdflist = []
    
    # get distribution of values in classes in a for loop
    for yearix in range(len(vals)):
        
        # get annual values
        yr = vals[yearix]
        
        # get uniques and occurrences
        #nuoc = np.unique(mc.Percentiles(yr, pct=[20,40,60,80,100]).yb, return_counts=True)
        nuoc = np.unique(mc.NaturalBreaks(yr, k=k, initial=40).yb, return_counts=True)
        
        # get top bins and lowest
        bins = np.around(mc.NaturalBreaks(yr, k=k, initial=40).bins, 3)
        lbin = round(yr.min(), 3)
        
        # produce dataframe
        classdf = pd.DataFrame({'class':nuoc[0],
                                'count':nuoc[1],
                                'bin_floor':lbin,
                                'bins':bins})
        
        # add year information
        classdf['year'] = yearix
        
        # append to list
        classdflist.append(classdf)
        
    # concatenate list
    cdf = pd.concat(classdflist, ignore_index=True)
    
     # plot
    fig, ax = plt.subplots(1,2, figsize=(13,5))
    ax = ax.flatten()
    g = sns.lineplot(data=cdf, x='year', y='count', hue='class',palette='muted',
                     ax=ax[0])
    ax[0].set_title('Class distribution with {}'.format(file))
    g = sns.lineplot(data=cdf, x='year', y='bins', hue='class', palette='muted',
                     ax=ax[1])
    ax[1].set_title('Class bin boundaries')
    plt.savefig(r'W:\maphel_langtime\plots\class_distribution_{}_natbre_k{}.pdf'.format(file,k),
                dpi=300, bbox_inches='tight')
    
    # get a markov chain object based on quintiles
    m5 = giddy.markov.Markov(q5)
    
    # get transition, probabiliy, average time for state change matrices for aspatial markov
    trans = m5.transitions
    probs = np.around(m5.p, decimals=5)
    avgs = giddy.ergodic.fmpt(m5.p) # mean first passage of time
    
    # plot spatial markov across spatial lags
    if k == 5:
        # set tick labels
        names = ['Very Low', 'Low', 'Moderate', 'High', 'Very High']

    elif k == 4:
        # set tick labels
        names = ['Low', 'Low Moderate', 'High Moderate', 'High']

    # get global probability matrix as pandas dataframe
    tdf = pd.DataFrame(m5.p, columns=names, index=names)
    tdf.to_pickle(r'W:\maphel_langtime\pickles\markov_probs_{}_natbre_k{}.pkl'.format(file,k))
    
    # get mean first passages as dataframe
    mfpt_df = pd.DataFrame(avgs, columns=names, index=names)
    mfpt_df.to_pickle(r'W:\maphel_langtime\pickles\markov_mfpt_{}_natbre_k{}.pkl'.format(file,k))
    
    # get spatial weight matrix for the grid KNN=8
    print('[INFO] - Calculating spatial Markov with KNN 8 weights matrix'\
          ' for {}...'.format(file))
    w = KNN.from_dataframe(df, k=8)
    
    # transpose values for spatial markov
    tvals = vals.transpose()
    
    # get spatial markov object
    #sm = giddy.markov.Spatial_Markov(tvals, w, fixed=True, k=k, m=k, variable_name='Quintiles')
    sm = giddy.markov.Spatial_Markov(q5, w, discrete=True,
                                     variable_name='Quintile',
                                     fill_empty_classes=True)
    
    # loop over spatial lag probabilities and save them
    for n, pm in enumerate(sm.P):
        
        #get spatial lag name and replace space with underscore
        sl = names[n].replace(' ','_')
        
        #convert to dataframe
        pmdf = pd.DataFrame(pm, columns=names, index=names)
        
        # save to pickle
        pmdf.to_pickle(r'W:\maphel_langtime\pickles\markov_SL_{}_natbre_probs_{}_k{}.pkl'.format(sl,file,k))
        
    
    # get statistical tests
    qstat = round(sm.Q, 3)
    qprob = "%.3f"%sm.Q_p_value
    lrstat = round(sm.LR, 3)
    lrprob = "%.3f"%sm.LR_p_value
    kullback = round(list(giddy.markov.kullback(sm.T).values())[0], 3)
    kbdof = list(giddy.markov.kullback(sm.T).values())[1]
    kbp = "%.3f"%list(giddy.markov.kullback(sm.T).values())[2]
    chi = round(giddy.markov.Homogeneity_Results(sm.T).Q, 3)
    cdof = giddy.markov.Homogeneity_Results(sm.T).dof
    cp = "%.3f"%giddy.markov.Homogeneity_Results(sm.T).Q_p_value
    
    # get a report text for stats
    report = 'Q-statistic: ' + str(qstat) + ' (p<' + qprob + ')                  '\
        + 'Likelihood-ratio: ' + str(lrstat) + ' (p<' + lrprob + ')\n'\
            'Kullback: ' + str(kullback) + ' (p<' + kbp + ' & DOF:' + str(kbdof) +')        ' \
                '\u03C7\u00b2: ' + str(chi) + ' (p<' + cp + ' & DOF:' + str(cdof) +')'
    
    # plot spatial markov across spatial lags
    if k == 5:
        fig, axes = plt.subplots(2,3, figsize=(16,9))
        # set tick labels
        ticks = ['Very Low', 'Low', 'Moderate', 'High', 'Very High']
        cnames = ['very low', 'low', 'moderate', 'high', 'very high']
        # flatten axes for easier looping
        axes = axes.flatten()
    elif k == 4:
        fig, axes = plt.subplots(2,3, figsize=(16,9))
        # set tick labels
        ticks = ['Low', 'Low\nModerate', 'High\nModerate', 'High']
        cnames = ['low', 'low moderate', 'high moderate', 'high']
        # flatten axes for easier looping
        axes = axes.flatten()
        # make last ax invisible
        axes[5].set_visible(False)

    
    # loop over axes
    for i, ax in enumerate(axes):
        # check if 5 classes
        if k == 5:
            if i == 0:
                # plot heatmap of global transition probabilities
                g = sns.heatmap(sm.p, annot=True, linewidths=.4, ax=ax, cbar=True,
                                vmin=0, vmax=1, square=True, cmap=cpal, fmt='.3f',
                                xticklabels=ticks, yticklabels=ticks)
                g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
                g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
                ax.set_title('Global', fontsize=12)
            else:
                # get transition probabilities for spatial lag
                p_temp = sm.P[i-1]
                # plot heatmap
                g = sns.heatmap(p_temp, annot=True, linewidths=.5, ax=ax, cbar=True, vmin=0,
                                vmax=1, square=True, cmap=cpal, fmt='.3f',
                                xticklabels=ticks, yticklabels=ticks)
                g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=9)
                g.set_xticklabels(g.get_xticklabels(), rotation=0, fontsize=9)
                ax.set_title('Surrounded by {} value class'.format(cnames[i - 1]), fontsize=11)
        # check if 4 classes
        elif k == 4:
            if i == 0:
                # plot heatmap of global transition probabilities
                g = sns.heatmap(sm.p, annot=True, linewidths=.4, ax=ax, cbar=True,
                                vmin=0, vmax=1, square=True, cmap=cpal, fmt='.3f',
                                xticklabels=ticks, yticklabels=ticks)
                g.set_yticklabels(g.get_yticklabels(), rotation=0)
                ax.set_title('Global', fontsize=12)
            elif i > 0 and i < 5:
                # get transition probabilities for spatial lag
                p_temp = sm.P[i-1]
                # plot heatmap
                g = sns.heatmap(p_temp, annot=True, linewidths=.5, ax=ax, cbar=True, vmin=0,
                                vmax=1, square=True, cmap=cpal, fmt='.3f',
                                xticklabels=ticks, yticklabels=ticks)
                g.set_yticklabels(g.get_yticklabels(), rotation=0)
                ax.set_title('Surrounded by {} value class'.format(cnames[i - 1]), fontsize=11)
            else:
                pass
    plt.figtext(0.34, 0.057, report, ha='left', fontsize=9)
    plt.savefig(r'W:\maphel_langtime\plots\spatial_markov_{}_natbre_k{}.pdf'.format(file,k),
                dpi=300, bbox_inches='tight')
    

