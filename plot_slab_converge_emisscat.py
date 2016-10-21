#!/usr/bin/env python
#
# Code to make a plot to illustrate the importance of scattering of the
#    dust emission
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 1 Aug 2016 - written
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

    # setup figure
    fig, ax = plt.subplots(ncols=2, figsize=(14,6))

    # factor for calibrating images into MJy/sr
    conv_fac = 159999.99
    im_conv_fac = 0.0039008563*conv_fac
    de_conv_fac = 0.39346823*conv_fac

    # get the stellar scattered component at 35.11 microns
    hdulist = fits.open('dirty_runs/slab_eff_grain/DIRTY_bm1_effgrain_tau_010.00_theta_090_w29_35.1119um.fits')
    stel_scat_image = hdulist['SCAT_SBoI'].data*im_conv_fac
    hdulist.close()

    # get the components of the dust emission at 35.11 microns
    hdulist = fits.open('dirty_runs/slab_eff_grain/DIRTY_bm1_effgrain_tau_010.00_theta_090_de_ge1_w29_35.1119um.fits')

    tot_image = hdulist['TOT_SBoI'].data*de_conv_fac
    scat_image = hdulist['SCAT_SBoI'].data*de_conv_fac
    direct_image = hdulist['STEL_SBoI'].data*de_conv_fac

    images = [tot_image+stel_scat_image, direct_image, scat_image, 
              stel_scat_image]
    labels = ['total','direct dust dust emission','scattered dust emission',
              'scattered stellar emission']

    hdulist.close()

    # info for the cut
    cut1 = np.array([95,115])

    syms = ['k-','k--','k-.','k:']

    # now process as before
    for i, timage in enumerate(images):

        # rot image
        timage = np.rot90(timage,1)

        timage = timage[19:520,90:510]
        if i == 0:
            nx = timage.shape[0]
            ny = timage.shape[1]
            cut1_plot_x = np.arange(nx)

        # first cut (y)
        cut1_plot_y = np.median(timage[:,cut1[0]:cut1[1]],axis=1)
        gindxs, = np.where(cut1_plot_y > 0.)
        if gindxs.shape[0] > 0:
            gindxs2 = gindxs[1:-1]

        ax[1].plot(cut1_plot_x[gindxs],cut1_plot_y[gindxs],
                     syms[i],label=labels[i])

    ax[1].set_xlim(80.,205)
    ax[1].set_yscale('log')

    #ax[1].set_ylabel('SB [MJy/sr]')
    ax[1].set_xlabel('pixel')

    ax[1].set_ylim(1e-9,1e2)

    ax[1].get_yaxis().set_ticks([])

    ax[1].set_title(r'DIRTY decomposed')

    leg = ax[1].legend(loc='upper left')
    leg.get_frame().set_linewidth(2)

    # now get and plot DIRTY, CRT, and TRADING
    dispnames = ['CRT','DART-ray','DIRTY','Hyperion','SKIRT',
                       'TRADING','SOC']
    modnames = ['crt','dartr','dirty','hyper','skirt','tradi','SOC']
    syms = ['b-','g-','r-','c-','m-','y-','k-','b--','g--','r--','c--',
               'm--','y--','k--']

    files = [modname + '/' + modname + 
             '_slab_eff_t1e1_i090a000_w035.11.fits'
             for modname in modnames]

    for i,cfile in enumerate(files):
        
        # read in the image
        timage = np.squeeze(fits.getdata(cfile))

        if dispnames[i] in ['CRT','SOC']:i
            timage = np.rot90(timage,3)

        if timage.shape[0] == 300:
            timage = timage[9:260,45:255]
        elif timage.shape[0] == 600:
            timage = timage[19:520,90:510]

        # first cut (y)
        cut1_plot_y = np.median(timage[:,cut1[0]:cut1[1]],axis=1)
        gindxs, = np.where(cut1_plot_y > 0.)
        if gindxs.shape[0] > 0:
            gindxs2 = gindxs[1:-1]

        ax[0].plot(cut1_plot_x[gindxs],cut1_plot_y[gindxs],
                   syms[i],label=dispnames[i])


    ax[0].set_xlim(80.,205)
    ax[0].set_yscale('log')

    ax[0].set_ylabel('SB [MJy/sr]')
    ax[0].set_xlabel('pixel')
    
    ax[0].set_ylim(1e-9,1e2)

    ax[0].set_title(r'Model Comparison')

    leg = ax[0].legend(loc='upper left')
    leg.get_frame().set_linewidth(2)

    # optimize the figure layout
    plt.tight_layout()

    # display the plot
    save_name = 'dirty_converge_emissscat'
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

