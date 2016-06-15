#!/usr/bin/env python
#
# Code to make a plot to show the different components of the SED
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 13 Jun 2016 - written
#

import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()

    parser.add_argument("--eps", help="Save the plot as an " + \
                        "encapsulated file",
                        action="store_true")
    parser.add_argument("--pdf", help="Save the plot as a " + \
                        "portable document file file",
                        action="store_true")
    parser.add_argument("--png", help="Save the plot as a " + \
                        "portable network graphics file",
                        action="store_true")
    args = parser.parse_args()

    # setup for nice plots
    fontsize = 24

    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.minor', width=2)

    # setup figure
    fig, ax = pyplot.subplots(figsize=(11,8.5))

    filename = 'dirty/dirty_slab_sto_t1e0_i150a000.sed'

    data = np.loadtxt(filename)
    # change the total dust emission column to be just direct dust
    #   emission instead of total dust emission
    data[:,4] -= data[:,5]

    # decmoposed sed components to display
    comp_indxs = [2,3,4,5,6,1]
    symtype = ['b-','b--','r-','r--','b:','k-']
    n_comps = len(comp_indxs)

    label_text = ['Direct Stellar','Scattered Stellar',
                  'Direct Dust Emission','Scattered Dust Emission',
                  'Transparent Stellar',
                  'Total']
    
    for k in range(n_comps):

        ck = comp_indxs[k]
        gindxs, = np.where(data[:,ck] > 0.0)
        ax.plot(data[gindxs,0],data[gindxs,ck],
                symtype[k],label=label_text[k])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    ax.set_ylabel('Flux [Jy]')
    ax.set_xlim([0.08,1.1e3])
    ax.set_ylim([1e-7,5e4])

    leg = ax.legend(loc=2,fontsize=0.8*fontsize, ncol=2)

    leg.get_frame().set_linewidth(2)

    # optimize the figure layout
    pyplot.tight_layout()

    # display the plot
    save_name = 'dirty_sed_comp_example'
    if args.png:
        fig.savefig(save_name+'.png')
        fig.savefig(save_name+'_small.png',dpi=11)
    elif args.eps:
        fig.savefig(save_name+'.eps')
    elif args.pdf:
        fig.savefig(save_name+'.pdf')
    else:
        pyplot.show()

    pyplot.close(fig)

