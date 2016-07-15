#!/usr/bin/env python
#
# Code to make a plot to make the nxy convergence figure
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 14 Jun 2016 - written
#

import argparse

import numpy as np
import matplotlib.pyplot as pyplot
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib

from astropy.io import fits

import cubehelix

from plot_slab_converge import plot_converge_slice
from plot_slab_converge_sed import plot_indiv_comp

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
    fontsize = 16

    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.minor', width=2)

    # setup figure
    fig, ax = pyplot.subplots(figsize=(12,8))

    # use gridspec to allow for one plot to be larger than the others
    gs = gridspec.GridSpec(2, 6, height_ratios=[2,1])

    ax = []
    ax.append(pyplot.subplot(gs[0,0:3]))
    ax.append(pyplot.subplot(gs[0,3:6]))
    ax.append(pyplot.subplot(gs[1,0]))
    ax.append(pyplot.subplot(gs[1,1]))
    ax.append(pyplot.subplot(gs[1,2]))
    ax.append(pyplot.subplot(gs[1,3]))
    ax.append(pyplot.subplot(gs[1,4]))
    ax.append(pyplot.subplot(gs[1,5]))

    plot_indiv_comp(ax[0], ['1e0','1e1'], ['000','180'], 'Total',
                    'dirty_nxy', r'$n_{xy}$', fontsize=fontsize)
    leg = ax[0].legend(fontsize=0.95*fontsize, loc='lower left')
    leg.get_frame().set_linewidth(2)
    ax[0].set_title('Global SED')
    ax[0].set_xlim(2,100)

    plot_converge_slice(ax[1], ['1e0','1e1'], ['035.11','151.99'], '000',
                        'dirty_nxy', r'$n_{xy}$', fontsize=fontsize)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    leg = ax[1].legend(fontsize=0.95*fontsize, loc='lower left')
    leg.get_frame().set_linewidth(2)
    ax[1].set_title(r'Y Image Slice, $\theta = 000^\circ$')
    ax[1].set_xlim(2,100)
    
    modname = 'dirty_nbinxy'
    cvals = ['5','10','20','50','100','200']

    custcmap = cubehelix.cmap(reverse = False, rot=1, start=0, 
                              startHue=320, sat = 2)

    minmax_y_vals = np.array([0.05,1])
    for i, cval in enumerate(cvals):
        tau = '1e0'
        angle = '000'
        wave = '035.11'
        cfile = modname + '/' + modname + \
            '_' + cval + \
            '_slab_eff_t' + tau + '_i' + \
            angle + 'a000_w' + \
            wave + '.fits'

        timage = np.squeeze(fits.getdata(cfile))

        timage = timage[45:255,45:255]

        cur_cax = ax[i+2].imshow(timage,
                                 norm=LogNorm(vmin=minmax_y_vals[0],
                                              vmax=minmax_y_vals[1]),
                                 origin='lower',
                                 cmap=custcmap)
        ax[i+2].get_xaxis().set_visible(False)
        ax[i+2].set_title(r'$n_{xy} = ' + cval + '$',
                        fontsize=fontsize)

        if i == 0:
            ax[i+2].set_ylabel(r'$\lambda = ' + wave + '$ $\mu m$')
            ax[i+2].yaxis.set_ticklabels([])
            ax[i+2].yaxis.set_ticks_position('none')
        else:
            ax[i+2].get_yaxis().set_visible(False)

            
    # optimize the figure layout
    gs.tight_layout(fig, w_pad=-0.25, h_pad=-0.25)

    # display the plot
    save_name = 'dirty_converge_nxy'
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

