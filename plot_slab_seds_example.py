#!/usr/bin/env python
#
# Code to make a plot to show the variation by observer angle theta
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 15 Jun 2016 - written
#

import argparse

import numpy as np
import matplotlib.pyplot as plt
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
    fig, ax = plt.subplots(figsize=(11,8.5))

    modname = 'dirty'

    taus = ['1e-2','1e-1','1e0','1e1']
    angles = ['000','090','180']

    coltype = ['b','g','m','r']
    symtype = ['-','--','-.',':']

    for j,tau in enumerate(taus):
        ifilenames = [modname + '/' + modname + '_slab_sto_t' + tau + '_i'+ angle + 'a000.sed'
                      for angle in angles]


        for k,cfile in enumerate(ifilenames):
            # read in the SED file for this model
            data = np.loadtxt(cfile)
    
            # change the total dust emission column to be just direct dust
            #   emission instead of total dust emission
            data[:,4] -= data[:,5]

            gindxs, = np.where(data[:,1] > 0.0)
            ax.plot(data[gindxs,0],data[gindxs,1],
                    coltype[j]+symtype[k],
                    label=r'$\tau_z = ' + tau + r';$ $\theta = ' + \
                        angles[k] + '$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    ax.set_ylabel('Flux [Jy]')
    ax.set_xlim([0.08,1.1e3])
    ax.set_ylim([1e-3,1.1e3])

    # Create two custom legends (more compact)

    # taus
    leg1 = ax.legend([plt.Line2D((0,1),(0,0), color='r', linestyle='-'),
                     plt.Line2D((0,1),(0,0), color='m', linestyle='-'),
                     plt.Line2D((0,1),(0,0), color='g', linestyle='-'),
                     plt.Line2D((0,1),(0,0), color='b', linestyle='-')],
                    [r'$\tau_z = 1e1$',
                     r'$\tau_z = 1e0$',
                     r'$\tau_z = 1e-1$',
                     r'$\tau_z = 1e-2$'],
                    fontsize=0.9*fontsize,
                    loc='upper left')
    leg1.get_frame().set_linewidth(2)

    # Add the legend manually to the current Axes.
    plt.gca().add_artist(leg1)

    # angles
    leg2 = ax.legend([plt.Line2D((0,1),(0,0), color='k', linestyle='-'),
                     plt.Line2D((0,1),(0,0), color='k', linestyle='--'),
                     plt.Line2D((0,1),(0,0), color='k', linestyle='-.')],
                    [r'$\theta = 0^\circ$',
                     r'$\theta = 90^\circ$',
                     r'$\theta = 180^\circ$'],
                    fontsize=0.9*fontsize,
                    loc='upper center')
    leg2.get_frame().set_linewidth(2)

    # optimize the figure layout
    plt.tight_layout()

    # display the plot
    save_name = 'dirty_sed_angles'
    if args.png:
        fig.savefig(save_name+'.png')
        fig.savefig(save_name+'_small.png',dpi=11)
    elif args.eps:
        fig.savefig(save_name+'.eps')
    elif args.pdf:
        fig.savefig(save_name+'.pdf')
    else:
        plt.show()

    plt.close(fig)

