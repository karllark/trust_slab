#!/usr/bin/env python
#
# Code to make a plot to make the nphot convergence figure
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 14 Jun 2016 - written
#

import argparse

import numpy as np
import matplotlib.pyplot as plt
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
    fontsize = 20

    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.minor', width=2)

    # setup figure
    fig, ax = plt.subplots(figsize=(18,12))

    # use gridspec to allow for one plot to be larger than the others
    gs = gridspec.GridSpec(4, 9, height_ratios=[0.435,0.130,0.435,1.0])

    ax = []
    ax.append(plt.subplot(gs[0:3,0:3]))
    ax.append(plt.subplot(gs[3,0:3]))
    ax.append(plt.subplot(gs[3,3:6]))
    ax.append(plt.subplot(gs[3,6:9]))

    ax.append(plt.subplot(gs[0,3]))
    ax.append(plt.subplot(gs[1,3]))
    ax.append(plt.subplot(gs[2,3]))

    ax.append(plt.subplot(gs[0,4]))
    ax.append(plt.subplot(gs[1,4]))
    ax.append(plt.subplot(gs[2,4]))

    ax.append(plt.subplot(gs[0,5]))
    ax.append(plt.subplot(gs[1,5]))
    ax.append(plt.subplot(gs[2,5]))

    ax.append(plt.subplot(gs[0,6]))
    ax.append(plt.subplot(gs[1,6]))
    ax.append(plt.subplot(gs[2,6]))

    ax.append(plt.subplot(gs[0,7]))
    ax.append(plt.subplot(gs[1,7]))
    ax.append(plt.subplot(gs[2,7]))

    ax.append(plt.subplot(gs[0,8]))
    ax.append(plt.subplot(gs[1,8]))
    ax.append(plt.subplot(gs[2,8]))

    plot_indiv_comp(ax[0], ['1e0','1e1'], ['000','090','180'], 'Total',
                    'dirty_nphot', r'$N$', fontsize=fontsize)

    # Create two custom legends (more compact)

    # taus
    leg1 = ax[0].legend([plt.Line2D((0,1),(0,0), color='k', linestyle='-'),
                         plt.Line2D((0,1),(0,0), color='k', linestyle='--')],
                        [r'$\tau_z = 1e0$',
                         r'$\tau_z = 1e1$'],
                        fontsize=fontsize,
                        loc='lower left')
    leg1.get_frame().set_linewidth(2)

    # Add the legend manually to the current Axes.
    ax[0].add_artist(leg1)

    # angles
    leg2 = ax[0].legend([plt.Line2D((0,1),(0,0), color='r', linestyle='-'),
                         plt.Line2D((0,1),(0,0), color='b', linestyle='-'),
                         plt.Line2D((0,1),(0,0), color='g', linestyle='-')],
                        [r'$\theta = 0^\circ$',
                         r'$\theta = 90^\circ$',
                         r'$\theta = 180^\circ$'],
                        fontsize=fontsize,
                        loc='lower right')
    leg2.get_frame().set_linewidth(2)

    #leg = ax[0].legend(fontsize=0.85*fontsize, loc='lower left',ncol=2)
    #leg.get_frame().set_linewidth(2)

    ax[0].set_title('Global SED')
    ax[0].set_ylim([1e-4,1e2])
    ax[0].set_xlabel('')
    ax[0].xaxis.set_ticklabels([])

    plot_converge_slice(ax[1], ['1e0','1e1'], 
                        ['000.15','000.53','035.11','151.99'], '000',
                        'dirty_nphot', r'$N$', fontsize=fontsize)
    #leg = ax[1].legend(fontsize=0.85*fontsize, loc='lower left', ncol=2)
    #leg.get_frame().set_linewidth(2)
    ax[1].set_title(r'Y Image Slice, $\theta = 000^\circ$')

    plot_converge_slice(ax[2], ['1e0','1e1'], 
                        ['000.15','000.53','035.11','151.99'], '090',
                        'dirty_nphot', r'$N$', fontsize=fontsize)

    #leg = ax[2].legend(fontsize=0.85*fontsize, loc='lower left', ncol=2)
    #leg.get_frame().set_linewidth(2)

    ax[2].yaxis.set_ticklabels([])
    ax[2].set_ylabel('')
    ax[2].set_title(r'Y Image Slice, $\theta = 090^\circ$')

    # taus
    leg1 = ax[3].legend([plt.Line2D((0,1),(0,0), color='k', linestyle='-'),
                         plt.Line2D((0,1),(0,0), color='k', linestyle='--')],
                        [r'$\tau_z = 1e0$',
                         r'$\tau_z = 1e1$'],
                        fontsize=fontsize,
                        loc='lower center', bbox_to_anchor=(-0.05,0.05))
    leg1.get_frame().set_linewidth(2)

    # Add the legend manually to the current Axes.
    ax[3].add_artist(leg1)

    # angles
    leg2 = ax[2].legend([plt.Line2D((0,1),(0,0), color='c', linestyle='-'),
                         plt.Line2D((0,1),(0,0), color='m', linestyle='-'),
                         plt.Line2D((0,1),(0,0), color='y', linestyle='-'),
                         plt.Line2D((0,1),(0,0), color='k', linestyle='-')],
                        [r'$\lambda = 000.15$',
                         r'$\lambda = 000.53$',
                         r'$\lambda = 035.11$',
                         r'$\lambda = 151.99$'],
                        fontsize=fontsize,
                        loc='lower center', bbox_to_anchor=(-0.05,0.05),
                        ncol=2)
    leg2.get_frame().set_linewidth(2)

    plot_converge_slice(ax[3], ['1e0','1e1'], 
                        ['000.15','000.53','035.11','151.99'], '180',
                        'dirty_nphot', r'$N$', fontsize=fontsize)

    #leg = ax[3].legend(fontsize=0.85*fontsize, loc='lower left', ncol=2)
    #leg.get_frame().set_linewidth(2)

    ax[3].yaxis.tick_right()
    ax[3].yaxis.set_label_position("right")
    ax[3].yaxis.set_ticks_position('both')
    ax[3].set_title(r'Y Image Slice, $\theta = 180^\circ$')

    custcmap = cubehelix.cmap(reverse = False, rot=1, start=0, 
                              startHue=320, sat = 2)

    modname = 'dirty_nphot'
    cvals = ['1e6','1e7','1e8']

    wave = '000.15'

    minmax_y_vals = np.array([[5e-4,5e-3],
                              [1e-45,1e0],
                              [1e-35,1e-30]])

    offval = 4
    for i, cval in enumerate(cvals):
        tau = '1e1'
        angles = ['000','090','180']

        for j, angle in enumerate(angles):

            cfile = modname + '/' + modname + \
                '_' + cval + \
                '_slab_eff_t' + tau + '_i' + \
                angle + 'a000_w' + \
                wave + '.fits'

            timage = np.squeeze(fits.getdata(cfile))
            
            if j == 1:
                timage = timage[90:240,90:510]
            else:
                timage = timage[90:510,90:510]

            axval = offval + i*3 + j
            cur_cax = ax[axval].imshow(timage,
                                       norm=LogNorm(vmin=minmax_y_vals[j,0],
                                                    vmax=minmax_y_vals[j,1]),
                                       origin='lower',
                                       cmap=custcmap)
            ax[axval].get_xaxis().set_visible(False)

            if j == 0:
                ax[axval].set_title(r'$N = ' + cval + '$',
                                    fontsize=fontsize)

            #if i == 0:
                #ax[axval].set_ylabel(r'$\theta = ' + angle + '^\circ$')
                #ax[axval].yaxis.set_ticklabels([])
                #ax[axval].yaxis.set_ticks_position('none')
            #else:
            ax[axval].get_yaxis().set_visible(False)

    wave = '035.11'

    minmax_y_vals = np.array([[1e-1,2],
                              [1e-7,1e1],
                              [1e-1,2]])

    offval = 4 + 3*3
    for i, cval in enumerate(cvals):
        tau = '1e1'
        angles = ['000','090','180']

        for j, angle in enumerate(angles):

            cfile = modname + '/' + modname + \
                '_' + cval + \
                '_slab_eff_t' + tau + '_i' + \
                angle + 'a000_w' + \
                wave + '.fits'

            timage = np.squeeze(fits.getdata(cfile))
            
            if j == 1:
                timage = timage[90:240,90:510]
            else:
                timage = timage[90:510,90:510]

            axval = offval + i*3 + j
            cur_cax = ax[axval].imshow(timage,
                                       norm=LogNorm(vmin=minmax_y_vals[j,0],
                                                    vmax=minmax_y_vals[j,1]),
                                       origin='lower',
                                       cmap=custcmap)
            ax[axval].get_xaxis().set_visible(False)

            if j == 0:
                ax[axval].set_title(r'$N = ' + cval + '$',
                                    fontsize=fontsize)

            if i == 2:
                ax[axval].set_ylabel(r'$\theta = ' + angle + '^\circ$')
                ax[axval].yaxis.set_label_position("right")
                ax[axval].yaxis.set_ticklabels([])
                ax[axval].yaxis.set_ticks_position('none')
            else:
                ax[axval].get_yaxis().set_visible(False)

            
    # add the labels for the wavelengths
    fig.text (0.5, 0.99, r'$\lambda = 000.15$ $\mu m$', 
              horizontalalignment='center',
              verticalalignment='top',fontsize=1.3*fontsize)

    fig.text (0.81, 0.99, r'$\lambda = 035.11$ $\mu m$', 
              horizontalalignment='center',
              verticalalignment='top',fontsize=1.3*fontsize)

    # optimize the figure layout
    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.98], w_pad=0.25, h_pad=0.5)

    # display the plot
    save_name = 'dirty_converge_nphot'
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

