#!/usr/bin/env python2.7
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

def plot_decompose_sed(modnames, moddisplaynames, tau, angle, save_eps=False, save_png=False):

    # generate the filename
    filenames = [modname + '/' + modname + '_slab_eff_t'+tau+'_i'
                 for modname in modnames]
    n_files = len(filenames)

    # plot information
    fig_label = r'Slab, $\tau (1 \mu m)$ = '+tau+r', $\theta$ = ' + angle
    symtype = ['b-','g-','r-','c-','m-','y-','k-']
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
        # filename
        cfile = cfile + angle + 'a000.sed'
        # if the '+' symbol was used (not requested, but...)
        if not os.path.isfile(cfile):
            cfile = cfile.replace('t1e','t1e+')

        # read in the SED file for this model
        data = np.loadtxt(cfile)

        # variable definitions based on the size of the 1st file read
        if i == 0:
            all_data = np.empty((data.shape[0],data.shape[1],n_files))
            ave_sed_comps = np.full((data.shape[0],data.shape[1]-1),0.0)
            ave_sed_comps_npts = np.full((data.shape[0],data.shape[1]-1),0.0)
            n_comps = data.shape[1] - 1

        # change the total dust emission column to be just direct dust emission instead of total dust emission
        data[:,4] -= data[:,5]
    
        # save the decomposed SEDs
        all_data[:,:,i] = data

        # sum for component
        for k in range(n_comps):
            gindxs, = np.where(data[:,comp_indxs[k]] > 0.0)
            if len(gindxs) > 0:
                ave_sed_comps[gindxs,k] += data[gindxs,comp_indxs[k]]
                ave_sed_comps_npts[gindxs,k] += 1.0

    # get the average by dividing by the number of models contributing to each wavelength point
    for k in range(n_comps):
        gindxs, = np.where(ave_sed_comps_npts[:,k] > 0)
        if len(gindxs) > 0:
            ave_sed_comps[gindxs,k] /= ave_sed_comps_npts[gindxs,k]

    # loop over the different components and plot the needed quantities
    for k in range(n_comps):

        # plot the SED and components to the big plot
        gindxs, = np.where(ave_sed_comps[:,k] > 0.0)
        if len(gindxs) > 0:
            ax[0].plot(all_data[gindxs,0,0],ave_sed_comps[gindxs,k],total_symtype[k],label=label_text[k])

        # plot the percentage difference for each model from the average
        for i in range(n_files):
            gindxs, = np.where(all_data[:,comp_indxs[k],i] > 0.0)
            if len(gindxs) > 0:
                y = 100.*(all_data[gindxs,comp_indxs[k],i] - ave_sed_comps[gindxs,k])/ave_sed_comps[gindxs,k]
                if comp_indxs[k] == 4:
                    # needed for to have a good legend
                    ax[k+1].plot(all_data[gindxs,0,i], y, symtype[i],label=moddisplaynames[i])
                else:
                    ax[k+1].plot(all_data[gindxs,0,i], y, symtype[i])

        # set the axis limits and type
        ax[k+1].set_xscale('log')
        ax[k+1].set_xlim([0.08,1.1e3])
        ax[k+1].set_ylabel('% difference')

        # label each % difference plot with the component name
        ylimits = ax[k+1].get_ylim()
        ax[k+1].text(7e2,ylimits[0]+0.85*(ylimits[1]-ylimits[0]),label_text[k],
                     fontsize=fontsize,horizontalalignment='right')

    # x axis label
    ax[3].set_xlabel('wavelength [$\mu$m]')
    ax[5].set_xlabel('wavelength [$\mu$m]')

    # label the big SED plot with the model details
    ylimits = ax[0].get_ylim()
    ax[0].text(1e-1,ylimits[0]+0.85*(ylimits[1]-ylimits[0]),fig_label,fontsize=1.3*fontsize)

    # big SED plot axis details
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Flux [Jy]')
    ax[0].set_xlim([0.08,1.1e3])
    ax[0].set_ylim([1e-10,5e2])

    # enable the two needed legends
    ax[0].legend(loc=3,fontsize=fontsize)
    ax[4].legend(loc=3,fontsize=fontsize)

    # optimize the figure layout
    pyplot.tight_layout()

    # display the plot
    save_name =  'slab_t' + tau + '_i' + angle + '_decomposed_sed_comp'

    if save_png:
        fig.savefig(save_name+'.png')
    elif save_eps:
        fig.savefig(save_name+'.eps')
    else:
        pyplot.show()

    pyplot.close(fig)


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tau", metavar='string', nargs=1,
                        help="Optical depth of model run [default=1e0]")
    parser.add_argument("-a", "--angle", metavar='string', nargs=1,
                        help="Viewing angle of model run [default=090]")
    parser.add_argument("-e", "--eps", help="Save the plot as an encapsulated file",
                        action="store_true")
    parser.add_argument("-p", "--png", help="Save the plot as a portable network graphics file",
                        action="store_true")
                        
    args = parser.parse_args()

    # use default or commandline option
    if args.angle:
        if 'all' in args.angle:
            angles = ['000','030','060','090','120','150','180']
        else:
            angles = args.angle
    else:
        angles = ['090']
        
    if args.tau:
        if 'all' in args.tau:
            taus = ['1e-2','1e-1','1e0','1e1']
        else:
            taus = args.tau[0]
    else:
        taus = ['1e0']

    # models to display
    moddisplaynames = ['DIRTY','SKIRT','SOC','TRADING']
    modnames = ['dirty','skirt','SOC','tradi']

    for angle in angles:
        for tau in taus:
            plot_decompose_sed(modnames, moddisplaynames, tau, angle, save_eps=args.eps, save_png=args.png)
