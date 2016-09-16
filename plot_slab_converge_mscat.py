#!/usr/bin/env python
#
# Code to make a plot to make the mscat convergence figure
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 16 Sep 2016 - written
#

import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib

from astropy.io import fits

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

    # control the maximum number of axis ticks
    max_yticks = 6
    yloc = plt.MaxNLocator(max_yticks)

    # setup figure
    fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(14,10), sharex=True)

    modname = 'dirty_mscat'
    cvals = ['1','2','5','10','500']
    pvals = cvals
    n_files = len(cvals)
    pcol = ['b','g','r','c','k']

    waves = ['000.15','035.11','151.99']

    for j, wave in enumerate(waves):
        cut1 = np.array([95,115])
        for i, cval in enumerate(cvals):
            tau = '1e0'
            angle = '090'
            cfile = modname + '/' + modname + \
                '_' + cval + \
                '_slab_eff_t' + tau + '_i' + \
                angle + 'a000_w' + \
                wave + '.fits'

            timage = np.squeeze(fits.getdata(cfile))

            if timage.shape[0] == 300:
                timage = timage[9:260,45:255]
            elif timage.shape[0] == 600:
                timage = timage[19:520,90:510]

            # save the image in a cube
            if i == 0:
                nx = timage.shape[0]
                ny = timage.shape[1]
                all_images = np.empty((timage.shape[0],timage.shape[1],n_files))
                minmax_vals = np.empty((n_files,2))
                if timage.shape[0] > 300:
                    cut1 *= 2

            all_images[:,:,i] = timage

        # get the cut to compare the rest to
        cut1_plot_x = np.arange(nx) * 2.0
        cut1_plot_y_all = np.median(all_images[:,cut1[0]:cut1[1],n_files-1],
                                    axis=1)

        for i, cval in enumerate(cvals):

            # first cut (y)
            cut1_plot_y = np.median(all_images[:,cut1[0]:cut1[1],i],axis=1)
            gindxs, = np.where(cut1_plot_y > 0.)
            if gindxs.shape[0] > 0:
                gindxs2 = gindxs[1:-1]

            ax[0,j].plot(cut1_plot_x[gindxs],cut1_plot_y[gindxs],pcol[i]+'-',
                         label=r'$m_\mathrm{scat} = ' + pvals[i] + '$')

            # percent difference plot for first cut
            y = 100.*((cut1_plot_y[gindxs] - cut1_plot_y_all[gindxs])/ 
                      cut1_plot_y_all[gindxs])
            ax[1,j].plot(cut1_plot_x[gindxs],y,pcol[i]+'-')

        ax[0,j].set_xlim(80.,205)
        ax[1,j].set_xlim(80.,205)

        ax[0,j].get_yaxis().set_ticks([])

        ax[0,j].set_yscale('log')

        if j == 0:
            ax[0,j].set_ylabel('SB [MJy/sr]')
            ax[1,j].set_ylabel('% difference')

        ax[1,j].set_ylim(-70.,10.)

        ax[1,j].yaxis.set_major_locator(yloc)

    # now tau = 10

    cvals = ['1','2','5','10','20','50','500']
    pvals = cvals
    n_files = len(cvals)

    for j, wave in enumerate(waves):
        cut1 = np.array([95,115])
        for i, cval in enumerate(cvals):
            tau = '1e1'
            angle = '090'
            cfile = modname + '/' + modname + \
                '_' + cval + \
                '_slab_eff_t' + tau + '_i' + \
                angle + 'a000_w' + \
                wave + '.fits'

            timage = np.squeeze(fits.getdata(cfile))

            if timage.shape[0] == 300:
                timage = timage[9:260,45:255]
            elif timage.shape[0] == 600:
                timage = timage[19:520,90:510]

            # save the image in a cube
            if i == 0:
                nx = timage.shape[0]
                ny = timage.shape[1]
                all_images = np.empty((timage.shape[0],timage.shape[1],n_files))
                minmax_vals = np.empty((n_files,2))
                if timage.shape[0] > 300:
                    cut1 *= 2

            all_images[:,:,i] = timage

        # get the cut to compare the rest to
        cut1_plot_x = np.arange(nx)
        cut1_plot_y_all = np.median(all_images[:,cut1[0]:cut1[1],n_files-1],
                                    axis=1)

        for i, cval in enumerate(cvals):

            # first cut (y)
            cut1_plot_y = np.median(all_images[:,cut1[0]:cut1[1],i],axis=1)
            gindxs, = np.where(cut1_plot_y > 0.)
            if gindxs.shape[0] > 0:
                gindxs2 = gindxs[1:-1]

            ax[2,j].plot(cut1_plot_x[gindxs],cut1_plot_y[gindxs],
                         label=r'$m_\mathrm{scat} = ' + pvals[i] + '$')

            # percent difference plot for first cut
            y = 100.*((cut1_plot_y[gindxs] - cut1_plot_y_all[gindxs])/ 
                      cut1_plot_y_all[gindxs])
            ax[3,j].plot(cut1_plot_x[gindxs],y)

        ax[2,j].set_xlim(80.,205)
        ax[3,j].set_xlim(80.,205)

        ax[2,j].get_yaxis().set_ticks([])

        ax[2,j].set_yscale('log')

        if j == 0:
            ax[2,j].set_ylabel('SB [MJy/sr]')
            ax[3,j].set_ylabel('% difference')

        ax[3,j].set_ylim(-100.,10.)

        ax[3,j].yaxis.set_major_locator(yloc)

    ax[3,0].set_xlabel('pixel')
    ax[3,1].set_xlabel('pixel')
    ax[3,2].set_xlabel('pixel')

    #ax[0,0].set_ylim(0.7e-6,1e2)
    #ax[0,1].set_ylim(1e-2,1e4)

    ax[0,0].set_title(r'$\lambda = 0.15$ $\mu m$')
    ax[0,1].set_title(r'$\lambda = 35.11$ $\mu m$')
    ax[0,2].set_title(r'$\lambda = 151.99$ $\mu m$')

    leg = ax[2,2].legend(loc='upper left',fontsize=fontsize*0.9, 
                         bbox_to_anchor=(0.02,1.6))
    leg.get_frame().set_linewidth(2)

    fig.text (0.97, 0.72, r'$\tau = 1; \theta = 90^\circ$', rotation=-90.0,
              horizontalalignment='center',
              verticalalignment='center',fontsize=1.5*fontsize)

    fig.text (0.97, 0.28, r'$\tau = 10; \theta = 90^\circ$', rotation=-90.0,
              horizontalalignment='center',
              verticalalignment='center',fontsize=1.5*fontsize)

    # optimize the figure layout
    plt.tight_layout(h_pad=0.0,w_pad=0.0,
                     rect=[0, 0, 0.97, 1.0])

    # display the plot
    save_name = 'dirty_converge_mscat'
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

