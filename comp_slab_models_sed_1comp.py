#!/usr/bin/env python
#
# Code to make plots for comparing the global decomposed SEDs
# for the TRUST BM1 Slab benchmark
#  code for a single component
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 22 Jan 2016 - written
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as pyplot

def plot_sed_1comp_core(ax, modnames, moddisplaynames, tau, angle, cindx, 
                        single_comp=-1, max_plot_diff=10.0):
    
    # generate the filename
    ifilenames = [modname + '_t' + tau + '_i'+ angle + 'a000.sed'
                 for modname in modnames]

    # check all the files exisit, adjust if not
    filenames = []
    displaynames = []
    fileindxs = []
    for i, cfile in enumerate(ifilenames):
        if os.path.isfile(cfile):
            filenames.append(cfile)
            displaynames.append(moddisplaynames[i])
            fileindxs.append(i)
        else:
            print(cfile + ' not found.')
    n_files = len(filenames)

    # plot information
    symtype = ['b-','g-','r-','c-','m-','y-','k-','b--','g--','r--','c--',
               'm--','y--','k--']
    total_symtype = ['k-','b--','b:','r--','r:','b--']
    fontsize = 12

    # decomposed sed components to display
    comp_indxs = [1,2,3,4,5,6]

    label_text = ['Total',
                  'Direct Stellar','Scattered Stellar',
                  'Direct Dust Emission','Scattered Dust Emission',
                  'Transparent Stellar']

    # read in the results from each model
    for i, cfile in enumerate(filenames):
        # read in the SED file for this model
        data = np.loadtxt(cfile)

        # variable definitions based on the size of the 1st file read
        if i == 0:
            all_data = np.empty((data.shape[0],data.shape[1],n_files))
            ave_sed_comps = np.full((data.shape[0],data.shape[1]-1),0.0)
            ave_sed_comps_npts = np.full((data.shape[0],data.shape[1]-1),0.0)
            n_comps = data.shape[1] - 1
            n_waves = data.shape[0]

        # change the total dust emission column to be just direct dust
        #   emission instead of total dust emission
        data[:,4] -= data[:,5]
    
        # save the decomposed SEDs
        all_data[:,:,i] = data

    # get the average by medianing the stack of results
    for k in range(n_comps):
        for j in range(n_waves):
            gindxs, = np.where(all_data[j,comp_indxs[k],:] > 0)
            if len(gindxs) >= 2:
                ave_sed_comps[j,k] = np.median(all_data[j,comp_indxs[k],gindxs])
            elif len(gindxs) >= 1:
                ave_sed_comps[j,k] = all_data[j,comp_indxs[k],gindxs[0]]

    # if single comp is set to a non-negative value, then use that model as
    #   the comparison instead of the average
    if single_comp >= 0:
        for k in range(n_comps):
            ave_sed_comps[:,k] = all_data[:,comp_indxs[k],single_comp]

    # loop over the different components and plot the needed quantities
    k = cindx

    # plot the percentage difference for each model from the average
    for i in range(n_files):
        gindxs, = np.where((all_data[:,comp_indxs[k],i] > 0.0) &
                           (ave_sed_comps[:,k] > 0.0))
        if len(gindxs) > 0:
            y = 100.*(all_data[gindxs,comp_indxs[k],i] -
                      ave_sed_comps[gindxs,k])/ave_sed_comps[gindxs,k]

        ax.plot(all_data[gindxs,0,i], y, symtype[fileindxs[i]])

        # set the axis limits and type
        ax.set_xscale('log')
        ax.set_xlim([0.08,1.1e3])
        ax.set_ylabel('% difference')

    # put a limit on the max % difference to show
    cur_ylim = ax.get_ylim()
    new_ylim = [max([cur_ylim[0],-1.0*max_plot_diff]),
                min([cur_ylim[1],max_plot_diff])]
    ax.set_ylim(new_ylim)

        
def plot_sed_1comp(modnames, moddisplaynames, cindx, 
                   single_comp=-1, max_plot_diff=10.0, save_str='',
                   save_eps=False, save_png=False, save_pdf=False,
                   plot_all=False):
    
    # angles and optical depths to display
    angles = ['000','030','060','090','120','150','180']
    taus = ['1e-2','1e-1','1e0','1e1']

    # plot information
    #fig_label = r'Slab, $\tau (1 \mu m)$ = '+tau+r', $\theta$ = ' + angle
    symtype = ['b-','g-','r-','c-','m-','y-','k-','b--','g--','r--','c--',
               'm--','y--','k--']
    total_symtype = ['k-','b--','b:','r--','r:','b--']
    fontsize = 12

    # decmoposed sed components to display
    comp_indxs = [1,2,3,4,5,6]

    label_text = ['Total',
                  'Direct Stellar','Scattered Stellar',
                  'Direct Dust Emission','Scattered Dust Emission',
                  'Transparent Stellar']

    # setup figure
    fig, ax = pyplot.subplots(nrows=len(angles), ncols=len(taus),
                              figsize=(20,15))

    # loop over the angles and taus and plot
    minmax_yaxis = np.empty((2),dtype=float)
    minmax_yaxis[0] = 1e30
    minmax_yaxis[1] = 1e-30
    for i, angle in enumerate(angles):
        for j, tau in enumerate(taus):
            plot_sed_1comp_core(ax[i,j], modnames, moddisplaynames, tau, angle,
                                cindx,
                                max_plot_diff=max_plot_diff)
            c_ylim = ax[i,j].get_ylim()
            if c_ylim[0] < minmax_yaxis[0]:
                minmax_yaxis[0] = c_ylim[0]
            if c_ylim[1] > minmax_yaxis[1]:
                minmax_yaxis[1] = c_ylim[1]

    for i, angle in enumerate(angles):
        for j, tau in enumerate(taus):
            ax[i,j].set_ylim(minmax_yaxis)

    # x axis label
    #ax[3].set_xlabel('wavelength [$\mu$m]')
    #ax[5].set_xlabel('wavelength [$\mu$m]')

    # optimize the figure layout
    pyplot.tight_layout()

    # display the plot
    save_name =  'slab_t' + tau + '_i' + angle + '_decomposed_sed_comp'
    
    if single_comp >= 0:
        save_name += '_scomp' + str(single_comp)

    if save_str != '':
        save_name += '_' + save_str

    if save_png:
        fig.savefig(save_name+'.png')
        fig.savefig(save_name+'_small.png',dpi=11)
    elif save_eps:
        fig.savefig(save_name+'.eps')
    elif save_pdf:
        fig.savefig(save_name+'.pdf')
    else:
        pyplot.show()

    pyplot.close(fig)

if __name__ == "__main__":

    print('Use comp_slab_models.py --sed_1comp')
