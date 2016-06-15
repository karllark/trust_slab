#!/usr/bin/env python
#
# Code to make a plot to show the different emission approximations
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

    gtypes = ['sto','equ','eff']
    #moddisplaynames = ['DI (gtype='+gtype+')' for gtype in gtypes]
    moddisplaynames = ['Full Solution','Equilibrium Only',
                       'Single Effective Grain']
    modnames = ['dirty','dirty','dirty']
    imodnames = ['dirty/dirty_slab_' + gtype for gtype in gtypes]

    tau = '1e0'
    angle = '090'
    ifilenames = [modname + '_t' + tau + '_i'+ angle + 'a000.sed'
                 for modname in imodnames]

    symtype = ['b-','g--','r-.','c-','m-','y-','k-','b--','g--','r--','c--',
               'm--','y--','k--']

    for k,cfile in enumerate(ifilenames):
        # read in the SED file for this model
        data = np.loadtxt(cfile)
    
        # change the total dust emission column to be just direct dust
        #   emission instead of total dust emission
        data[:,4] -= data[:,5]

        gindxs, = np.where(data[:,1] > 0.0)
        ax.plot(data[gindxs,0],data[gindxs,1],
                symtype[k],label=moddisplaynames[k])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    ax.set_ylabel('Flux [Jy]')
    ax.set_xlim([0.08,1.1e3])

    leg = ax.legend(loc=2,fontsize=0.95*fontsize)

    leg.get_frame().set_linewidth(2)

    # optimize the figure layout
    pyplot.tight_layout()

    # display the plot
    save_name = 'dirty_emis_type'
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

