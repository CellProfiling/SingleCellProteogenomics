#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 10:41:42 2018

@author: devinsullivan
"""

import numpy as np
import matplotlib.pyplot as plt
import operator
import os
import seaborn as sns
import pandas as pd
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable

#from continuous_time_polar import psin_predict
def psin_predict(p,x):
    return p[0]*np.sin(np.pi*x**p[1])**p[2]+p[3]

def temporal_heatmap(model_free_list, plate_well_to_ab, loc_set,outfolder):
    #PROTEINS SELECTED FOR HIGHLIGHTING
    HIGHLIGHTS = ['ORC6','CCNE1','PHLDB1','DPH2']
    HIGHLIGHTS_MINOR = [ 'MCM10', 'ZNF32', 'JUN',  
                  'DUSP19', 'CCNB1', 'AURKB', 'BUB1B', 'PAPSS1',
                  'N6AMT1', 'FLI1']
    #HIGHLIGHTS['ORC6','MCM10','ZNF32','JUN','CCNE1','DUSP19',
    #   'CCNB1','AURKB','BUB1B','PAPSS1','N6AMT1','PHLDB1','DPH2','FLI1']
    #TIMING OF PHASE TRANSITIONS (MANUALLY DETERMINED BY DIANA)
    #hours (for the G1/S cutoff)
    G1_LEN = 10.833
    #hours (plus 10.833, so 13.458hrs for the S/G2 cutoff)
    G1_S_TRANS = 2.625
    #hours (plus 10.833 and 2.625 so 25.433 hrs for the G2/M cutoff)
    S_G2_LEN = 11.975
    #hours (this should be from the G2/M cutoff above to the end)
    #M_LEN = 0.5
    #We are excluding Mphase from this analysis
    TOT_LEN = G1_LEN+G1_S_TRANS+S_G2_LEN
    G1_PROP = G1_LEN/TOT_LEN
    G1_S_PROP = G1_S_TRANS/TOT_LEN+G1_PROP
    S_G2_PROP = S_G2_LEN/TOT_LEN+G1_S_PROP

    #create localization matrix
    #sort the coefficients based on the max value point
    loc_mat = np.zeros([len(model_free_list),len(loc_set)])
    for i,response in enumerate(model_free_list):
        curr_well = response[0]
        print(curr_well)
        print(plate_well_to_ab[curr_well])
        curr_locs = plate_well_to_ab[curr_well][3]
        loc_mat[i,curr_locs] = 1


    model_free_list.sort(key=operator.itemgetter(3))
    prob_interact = np.zeros([len(model_free_list),len(model_free_list)])
    sum_rows = []
    prot_list = []
    sorted_gene_array = []
    for i,response in enumerate(model_free_list):
        sorted_gene_array.append(response[4])
        curr_well = response[0]
        #set all other rows that share a loc to 1
        curr_locs = loc_mat[i,:]
        match_rows = np.greater(np.sum(loc_mat*curr_locs,axis=1),0)*response[3]
        sum_rows.append(np.sum(match_rows))
        prob_interact[i,:] = match_rows
        prot_list.append(response[1])

    fig, ax = plt.subplots(figsize=(10, 10))
    sc = ax.imshow(sorted_gene_array, interpolation='nearest')
    

    #Arange the x ticks
    xtick_labels = [str(np.around(x,decimals=2)) for x in np.linspace(0,1,11)]
    xtick_labels = xtick_labels#+['G1/S','S/G2']
    xphase_labels = ['G1/S','S/G2']
    my_xticks = np.arange(-.5, 20, 2)
    num_ticks = 20

    #my_xticks = np.concatenate(
    #        [my_xticks,
    #         np.asarray([G1_PROP*num_ticks-0.5, G1_S_PROP*num_ticks-0.5])])
    phase_trans = np.asarray([G1_PROP*num_ticks-0.5, G1_S_PROP*num_ticks-0.5])

    ax.set_xticks(my_xticks,minor=True)
    ax.set_xticklabels(xtick_labels,minor=True)
    ax.set_xticks(phase_trans, minor=False)
    ax.set_xticklabels(xphase_labels, minor=False)
    ax.tick_params(length=12)
    plt.xticks(size=12,fontname='Arial')
    plt.yticks(size=10,fontname='Arial')
#    ax.set_xticks([G1_PROP*num_ticks,G1_S_PROP*num_ticks,S_G2_PROP*num_ticks],minor=False)
#    ax.set_xticklabels(['g1','g1s','sg2'])
#

    #Do the y ticks
    ytick_locs = [prot_list.index(p) for p in HIGHLIGHTS]
    ytick_locs_minor = [prot_list.index(p) for p in HIGHLIGHTS_MINOR]
    ax.set_yticks(ytick_locs,minor=False)
    ax.set_yticklabels(HIGHLIGHTS,minor=False)
    ax.set_yticks(ytick_locs_minor,minor=True)
    ax.set_yticklabels(HIGHLIGHTS_MINOR,minor=True)
    ax.tick_params(direction='out', length=12, width=2, colors='k', 
                   axis='x',which='major')
    ax.tick_params(direction='out', length=56, width=1, colors='k', 
                   axis='y',which='major')

    ax.set_aspect('auto')
    plt.xlabel('Normalized division cycle',size=20,fontname='Arial')
    plt.ylabel('Gene',size=20,fontname='Arial')
    
    divider1 = make_axes_locatable(ax)
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(sc,cax = cax1)
    cbar.set_label('Relative expression',fontname='Arial',size=20)
    cbar.ax.tick_params(labelsize=18)
    
    plt.show()
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'sorted_heatmap21_sw30_take4.pdf'), transparent=True)
    return fig,ax,prot_list

def p_v_fracvar(plist_tot,perc_var_tot,fdr_tot,outfolder):
#DO A BUNCH OF PLOTTING
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.scatter(np.log10(plist_tot),perc_var_tot,c=fdr_tot)
    plt.xlabel('log10 pval',size=20,fontname='Arial')
    plt.ylabel('frac explained var',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'p_vs_frac_var_explained_compartments.pdf'))

def gini_v_fracvar(gini_tot,perc_var_tot,fdr_tot,outfolder):
    ##MT CONSTANTS -- Calculated by running in mt mode, here we use cytosol
    MT_GINI = 0.06600124369521103
    MT_PERC_VAR = 0.09086819*100
    #MT_FDR = 0.3783993#FROM DIANA: cyto 0.3783993  and cell 0.0756249

    #if it's a fraction, convert it to a percentage
    if np.max(perc_var_tot)<=1:
        perc_var_tot = perc_var_tot*100

    plt.close()
    fig, ax = plt.subplots(figsize=(10, 10))
    n_log_fdr = np.log10(fdr_tot)
    plt.scatter(gini_tot[::-1],perc_var_tot[::-1],c=n_log_fdr[::-1])
    plt.plot([np.min(gini_tot),np.max(gini_tot)],
              [np.mean(perc_var_tot),np.mean(perc_var_tot)],c='k')
    plt.plot([np.min(gini_tot),np.max(gini_tot)],
              [np.mean(perc_var_tot)+np.std(perc_var_tot),
               np.mean(perc_var_tot)+np.std(perc_var_tot)],c='k',linestyle='dashed')
    plt.scatter(MT_GINI,
              MT_PERC_VAR,c='r',s=100)
    plt.xlabel('gini',size=20,fontname='Arial')
    plt.ylabel('% explained variance',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'gini_vs_frac_var_explained_compartments.pdf'))

def gini_v_fracvar_2(all_coeff_list,outfolder):
    #SET THE INDS FOR DIFFERENT PARAMETERS (from continuous_time_polar.py)
    #Should really make this an object...would solve a lot of issues...

    ##MT CONSTANTS -- Calculated by running in mt mode, here we use cytosol
    MT_GINI = 0.06600124369521103
    MT_PERC_VAR = 0.09086819*100
    #MT_FDR = 0.3783993#FROM DIANA: cyto 0.3783993  and cell 0.0756249
    #PROTEINS SELECTED FOR HIGHLIGHTING
    HIGHLIGHTS = ['DUSP18','CCNB1','ORC6','KIF23','FAM71F1']
    #AXIS OFFSETS FOR ANNOTATIONS
    XOFFSET = 0.005
    YOFFSET = -0.005

    #make the lists
    perc_var_tot = []
    gini_tot = []
    fdr_tot = []
    H_inds = []
    H_labels = []
    for i,item in enumerate(all_coeff_list):
        protname = item.outsuff.split('_')[0]
        pv = item.perc_var_compartment
        gini = item.gini_compartment
        fdr = item.curr_fdr
        if protname in HIGHLIGHTS:
            H_inds.append(i)
            H_labels.append(protname)
        perc_var_tot.append(pv)
        gini_tot.append(gini)
        fdr_tot.append(fdr)

    #if it's a fraction, turn it to a percentage
    if np.max(perc_var_tot)<=1:
        perc_var_tot = np.asarray(perc_var_tot)*100


    #do the plots
    plt.close()
    fig, ax = plt.subplots(figsize=(10, 10))
    n_log_fdr = np.log10(fdr_tot)
    sc = ax.scatter(gini_tot[::-1],perc_var_tot[::-1],c=n_log_fdr[::-1])
    cbar = plt.colorbar(sc)
    cbar.set_label(r'$\ -log_{10}(FDR)$',fontname='Arial',size=24)
    cbar.ax.tick_params(labelsize=18)

    for ind,txt in zip(H_inds,H_labels):
        ax.annotate(txt,(gini_tot[ind]+XOFFSET,perc_var_tot[ind]+YOFFSET))

    ax.plot([np.min(gini_tot),np.max(gini_tot)],
              [np.mean(perc_var_tot),np.mean(perc_var_tot)],c='k',
              label=r'$\ \%var_{\mu}='+str(np.round(
                      np.mean(perc_var_tot),decimals=0)))
    ax.scatter(MT_GINI,
              MT_PERC_VAR,c='r',s=100,label='Microtubules')
    ax.annotate('MT',(MT_GINI+XOFFSET,MT_PERC_VAR+YOFFSET))
    plt.xlabel('Gini score',size=24,fontname='Arial')
    plt.ylabel('% Explained variance',size=24,fontname='Arial')
    plt.xticks(size=18,fontname='Arial')
    plt.yticks(size=18,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'gini_vs_frac_var_explained_compartments_2.pdf'))

def var_v_fracvar(var_tot,perc_var_tot,fdr_tot,outfolder):
    plt.close()
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.scatter(var_tot,perc_var_tot,c=fdr_tot)
    plt.xlabel('var',size=20,fontname='Arial')
    plt.ylabel('frac explained var',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'var_vs_frac_var_explained_compartments.pdf'))

def fracvar_compartments(perc_varlist,nc_perc_varlist,outfolder):
    plt.close()
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.plot(perc_varlist,c='b')
    plt.plot(nc_perc_varlist,c='r')
    plt.xlabel('gene',size=20,fontname='Arial')
    plt.ylabel('frac explained var',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'frac_var_explained_compartments.pdf'))

def fdr_v_fracvar(fdr_tot,perc_var_tot,outfolder):
    fig, ax = plt.subplots(figsize=(10, 10))
    n_log_fdr = -np.log10(fdr_tot)
    if np.max(perc_var_tot)<=1:
        perc_var_tot = perc_var_tot*100
    plt.scatter(n_log_fdr,perc_var_tot,c='grey')
    plt.xlabel(r'$\ -log_{10}(FDR)$',size=20,fontname='Arial')
    plt.ylabel('% explained variance',size=20,fontname='Arial')
    plt.plot([np.min(n_log_fdr),np.max(n_log_fdr)],
              [np.mean(perc_var_tot),np.mean(perc_var_tot)],c='r',
              label=r'$\ \%var_{\mu}=$'+str(int(np.round(np.mean(perc_var_tot),decimals=0)))+'%')
    plt.plot([np.min(n_log_fdr),np.max(n_log_fdr)],
              [10,
               10],c='k',label='10% var')
    plt.plot([-np.log10(0.05),-np.log10(0.05)],
              [np.min(perc_var_tot),np.max(perc_var_tot)],c='c',
              label='FDR=0.05')
    plt.xticks(size=18,fontname='Arial')
    plt.yticks(size=18,fontname='Arial')
    plt.legend(loc=2,prop={'size': 20})
    plt.tight_layout()
    plt.savefig(outfolder+os.sep+'fdr_vs_frac_var.pdf')

def temporal_mov_avg(curr_pol,curr_ab_norm,x_fit,y_fit,outname,outsuff):
    plt.close()
    #plot data
    plt.scatter(curr_pol,curr_ab_norm,c='c')
    plt.plot(x_fit,y_fit,'m')
    #label stuff
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.xlabel('pseudo-time')
    plt.ylabel(outsuff+' expression')
    plt.tight_layout()
    outfile = os.path.join(outname,outsuff+'_mvavg.pdf')
    if not os.path.exists(os.path.dirname(outfile)):
        os.path.makedirs(os.path.dirname(outfile))
    plt.savefig(outfile)

def temporal_psin(curr_pol,curr_ab_norm,fit_coeffs,outname,outsuff,n_pts=200):
    x_fit = np.linspace(0,1,n_pts)
    plt.close()
    fig, ax = plt.subplots(figsize=(10, 10),facecolor='None')
    #plot data
    plt.scatter(curr_pol,curr_ab_norm,c='grey')
    #plot fit
    y_fit = psin_predict(fit_coeffs.x,x_fit)
    plt.plot(x_fit,y_fit,'m')
    #label stuff
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.xticks(size=12,fontname='Arial')
    plt.yticks(size=12,fontname='Arial')
    plt.xlabel('pseudo-time',size=20,fontname='Arial')

    plt.ylabel(str.split(outsuff,'_')[0]+' normalized expression', size=20, fontname='Arial')
    plt.tight_layout()
    outfile = os.path.join(outname,outsuff+'_psin.pdf')
    if not os.path.exists(os.path.dirname(outfile)):
        os.path.makedirs(os.path.dirname(outfile))
    plt.savefig(outfile)

def hist2d_time(green_fucci,red_fucci,
                curr_fgreen,curr_fred,
                curr_pol,outname,outsuff,nbins=150):
    my_cmap = plt.cm.gray_r
    my_cmap.set_under(color='w',alpha=None)
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect(1)
    ax.hist2d(np.log10(green_fucci),
               np.log10(red_fucci),bins=nbins,cmap=my_cmap)
    cfgreen_log = np.log10(curr_fgreen)
    cfred_log = np.log10(curr_fred)
    #plt.plot(cfgreen_log,cfred_log,c='r',linewidth=0.5)
    #divider1 = make_axes_locatable(ax)
    #cax1 = divider1.append_axes("right", size="5%", pad=0.05)

    sc = ax.scatter(cfgreen_log,cfred_log,c=np.round(curr_pol,1))

    #cbar = plt.colorbar(sc,cax = cax1)
#    cbar = plt.colorbar(sc)
#    cbar.set_label('pseudo-time',fontname='Arial',size=24)
#    cbar.ax.tick_params(labelsize=18)
    #For internal use
    #label with the percent variance accounted for in the fucci channels
    #plt.xlabel('fgreen: '+str(np.round(perc_var_fgreen,2)))
    #plt.ylabel('fred: '+str(np.round(perc_var_fred,2)))

    #plt.title(outsuff)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_label('pseudo-time',fontname='Arial',size=24)
    cbar.ax.tick_params(labelsize=18)

    plt.sca(ax)
    plt.xlabel(r'$\propto log_{10}(GMNN_{fucci})$',size=20,fontname='Arial')
    plt.ylabel(r'$\propto log_{10}(CDT1_{fucci})$',size=20,fontname='Arial')
    plt.xticks(size=18,fontname='Arial')
    plt.yticks(size=18,fontname='Arial')
    #fig.tight_layout()
    plt.savefig(os.path.join(outname,outsuff+'_fucci2.pdf'))

def peak_time_v_perc_var(all_coeff_list,outfolder):
    #PROTEINS SELECTED FOR HIGHLIGHTING
    HIGHLIGHTS = ['DUSP18','CCNB1','ORC6','KIF23','FAM71F1']
    #AXIS OFFSETS FOR ANNOTATIONS
    XOFFSET = 0.005
    YOFFSET = -0.005
    #TIMING OF PHASE TRANSITIONS (MANUALLY DETERMINED BY DIANA)
    #hours (for the G1/S cutoff)
    G1_LEN = 10.833
    #hours (plus 10.833, so 13.458hrs for the S/G2 cutoff)
    G1_S_TRANS = 2.625
    #hours (plus 10.833 and 2.625 so 25.433 hrs for the G2/M cutoff)
    S_G2_LEN = 11.975
    #hours (this should be from the G2/M cutoff above to the end)
    #M_LEN = 0.5
    #We are excluding Mphase from this analysis
    TOT_LEN = G1_LEN+G1_S_TRANS+S_G2_LEN
    G1_PROP = G1_LEN/TOT_LEN
    G1_S_PROP = G1_S_TRANS/TOT_LEN+G1_PROP
    S_G2_PROP = S_G2_LEN/TOT_LEN+G1_S_PROP

    #make the lists
    perc_var_tot = []
    peaks_tot = []
    H_inds = []
    H_labels = []
    for i,item in enumerate(all_coeff_list):
        protname = item.outsuff.split('_')[0]
        pv = item.perc_var_compartment
        peak_loc = item.max_val_compartment
        if protname in HIGHLIGHTS:
            H_inds.append(i)
            H_labels.append(protname)
        perc_var_tot.append(pv)
        peaks_tot.append(peak_loc)

    #do the plots
    plt.close()
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(peaks_tot,perc_var_tot)
    for ind,txt in zip(H_inds,H_labels):
        ax.annotate(txt,(peaks_tot[ind]+XOFFSET,perc_var_tot[ind]+YOFFSET))
    ax.plot([np.min(peaks_tot),np.max(peaks_tot)],
              [np.mean(perc_var_tot),np.mean(perc_var_tot)],c='k')
    ax.plot([np.min(peaks_tot),np.max(peaks_tot)],
              [np.mean(perc_var_tot)+np.std(perc_var_tot),
               np.mean(perc_var_tot)+np.std(perc_var_tot)],c='k',linestyle='dashed')

    plt.xlabel('peak_pos',size=20,fontname='Arial')
    numeric_ticks = np.round(np.linspace(0,1,11),decimals=1)
    ax.set_xticks(numeric_ticks,minor=True)
    ax.set_xticklabels(numeric_ticks,minor=True)
    xphase_labels = ['G1/S','S/G2','G2/M']
    phase_trans = [G1_PROP,G1_S_PROP,S_G2_PROP]
    ax.set_xticks(phase_trans, minor=False)
    ax.set_xticklabels(xphase_labels, minor=False)

    plt.ylabel('frac explained var',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'peakpos_vs_frac_var_explained_compartments.pdf'))

    plt.close()
    fig, ax = plt.subplots(figsize=(10, 10))
    prev_phase = 0
    phase_data = np.zeros(np.shape(peaks_tot))
    var_data = []
    for i, (phase, plabel) in enumerate(zip(phase_trans,xphase_labels)):
        curr_inds = [np.greater_equal(peaks_tot,prev_phase) * np.less(peaks_tot,phase)]

        phase_data[curr_inds] = i
#        var_data.append([perc_var_tot[curr_ind]])

        prev_phase = phase

    phase_df = pd.DataFrame()
    phase_df = phase_df.assign(phase = phase_data)
    phase_df = phase_df.assign(perc_var = perc_var_tot)

    g1_data = phase_df.loc[phase_df['phase'] == 0]
    s_data = phase_df.loc[phase_df['phase'] == 1]
    g2_data = phase_df.loc[phase_df['phase'] == 2]

    [tstat_g1s,pval_g1s] = stats.ttest_ind(g1_data['perc_var'],s_data['perc_var'])
    [tstat_sg2,pval_sg2] = stats.ttest_ind(s_data['perc_var'],g2_data['perc_var'])
    [tstat_g1g1,pval_g1g2] = stats.ttest_ind(g1_data['perc_var'],g2_data['perc_var'])


    sns.violinplot(x='phase',y='perc_var',data=phase_df)

def fucci_hist2d(centered_data,cart_data_ur,start_pt,outfolder,nbins=200):
    fig, ax1 = plt.subplots(figsize=(10,10))
    mycmap = plt.cm.gray_r
    mycmap.set_under(color='w',alpha=None)
    ax1.hist2d(centered_data[:,0],centered_data[:,1],bins=nbins,alpha=1,cmap=mycmap)
    hist, xbins, ybins = np.histogram2d(cart_data_ur[0],cart_data_ur[1], bins=nbins, normed=True)
    extent = [xbins.min(),xbins.max(),ybins.min(),ybins.max()]
    im = ax1.imshow(
            np.ma.masked_where(hist == 0, hist).T,
            interpolation='nearest',
            origin='lower',
            extent=extent,
            cmap='plasma')
    plt.scatter(start_pt[0],start_pt[1],c='c',linewidths=4)
    plt.scatter(0,0,c='m',linewidths=4)
    plt.xlabel(r'$\propto log_{10}(GMNN_{fucci})$',size=20,fontname='Arial')
    plt.ylabel(r'$\propto log_{10}(CDT1_{fucci})$',size=20,fontname='Arial')
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,'masked_polar_hist.pdf'),transparent=True)
