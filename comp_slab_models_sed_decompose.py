#!/usr/bin/env python
#
# Code to make plots for comparing the global decomposed SEDs
# for the TRUST BM1 Slab benchmark
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 1 Jul 2015 - written
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec

def plot_decompose_sed(modnames, moddisplaynames, tau, angle,
                       single_comp=-1, max_plot_diff=10.0, save_str='',
                       save_eps=False, save_png=False, save_pdf=False,
                       plot_all=False):

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
    fig_label = r'Slab, $\tau (1 \mu m)$ = '+tau+r', $\theta$ = ' + angle
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
    fig, ax = pyplot.subplots(figsize=(15,10))

    # use gridspec to allow for one plot to be larger than the others
    gs = gridspec.GridSpec(4, 2)
    ax = [pyplot.subplot(gs[0:2,0]),
          pyplot.subplot(gs[0,1]),
          pyplot.subplot(gs[2,0]),
          pyplot.subplot(gs[3,0]),
          pyplot.subplot(gs[2,1]),
          pyplot.subplot(gs[3,1]),
          pyplot.subplot(gs[1,1])]
      
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
                #print(all_data[j,k,gindxs])
                #print(ave_sed_comps[j,k])
            elif len(gindxs) >= 1:
                ave_sed_comps[j,k] = all_data[j,comp_indxs[k],gindxs[0]]

    # if single comp is set to a non-negative value, then use that model as
    #   the comparison instead of the average
    if single_comp >= 0:
        for k in range(n_comps):
            ave_sed_comps[:,k] = all_data[:,comp_indxs[k],single_comp]

    # loop over the different components and plot the needed quantities
    for k in range(n_comps):

        # plot the SED and components to the big plot
        gindxs, = np.where(ave_sed_comps[:,k] > 0.0)
        if len(gindxs) > 0:
            ax[0].plot(all_data[gindxs,0,0],ave_sed_comps[gindxs,k],
                       total_symtype[k],label=label_text[k])
            if plot_all:
                for z in range(n_files):
                    ax[0].plot(all_data[gindxs,0,0],
                               all_data[gindxs,comp_indxs[k],z],total_symtype[k])

        # plot the percentage difference for each model from the average
        for i in range(n_files):
            gindxs, = np.where((all_data[:,comp_indxs[k],i] > 0.0) &
                               (ave_sed_comps[:,k] > 0.0))
            if len(gindxs) > 0:
                y = 100.*(all_data[gindxs,comp_indxs[k],i] -
                          ave_sed_comps[gindxs,k])/ave_sed_comps[gindxs,k]
                if (comp_indxs[k] == 4) & (i < 6):
                    # needed for to have a good legend
                    ax[k+1].plot(all_data[gindxs,0,i], y, symtype[fileindxs[i]],
                                 label=displaynames[i])
                if (comp_indxs[k] == 5) & (i >= 6):
                    ax[k+1].plot(all_data[gindxs,0,i], y, symtype[fileindxs[i]],
                                 label=displaynames[i])
                else:
                    ax[k+1].plot(all_data[gindxs,0,i], y, symtype[fileindxs[i]])

        # set the axis limits and type
        ax[k+1].set_xscale('log')
        ax[k+1].set_xlim([0.08,1.1e3])
        ax[k+1].set_ylabel('% difference')

        # put a limit on the max % difference to show
        cur_ylim = ax[k+1].get_ylim()
        new_ylim = [max([cur_ylim[0],-1.0*max_plot_diff]),
                    min([cur_ylim[1],max_plot_diff])]
        ax[k+1].set_ylim(new_ylim)

        # label each % difference plot with the component name
        ylimits = ax[k+1].get_ylim()
        ax[k+1].text(7e2,ylimits[0]+0.85*(ylimits[1]-ylimits[0]),label_text[k],
                     fontsize=fontsize,horizontalalignment='right')

    # x axis label
    ax[3].set_xlabel('wavelength [$\mu$m]')
    ax[5].set_xlabel('wavelength [$\mu$m]')

    # big SED plot axis details
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Flux [Jy]')
    ax[0].set_xlim([0.08,1.1e3])
    cur_ylim = ax[0].get_ylim()
    new_ylim = [cur_ylim[1]/1e11, cur_ylim[1]]
    ax[0].set_ylim(new_ylim)

    # label the big SED plot with the model details
    ylimits = ax[0].get_ylim()
    ypos = 10**(np.log10(ylimits[0])+0.92*(np.log10(ylimits[1])-
                                           np.log10(ylimits[0])))
    ax[0].text(1e-1,ypos,fig_label,fontsize=1.3*fontsize)

    # enable the two needed legends
    ax[0].legend(loc=3,fontsize=fontsize)
    ax[4].legend(loc=2,fontsize=fontsize)
    if n_files > 6:
        ax[5].legend(loc=2,fontsize=fontsize)

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

    print('Use comp_slab_models.py')
