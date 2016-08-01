#!/usr/bin/env python
#
# Code to make a plot to show the different images for the Slab model
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 13 Jun 2016 - written
#

import argparse

import numpy as np
import matplotlib.pyplot as pyplot
from matplotlib.colors import LogNorm
import matplotlib

from astropy.io import fits

import cubehelix

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()

    good_waves = ['000.15','035.11','151.99']
    parser.add_argument("-w", "--wave", choices=good_waves, default='000.15',
                        help="wavelength to display")

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
    fontsize = 18

    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.minor', width=2)

    # setup figure
    fig, ax = pyplot.subplots(nrows=3,ncols=5,figsize=(10,6), 
                              gridspec_kw = 
                              {'height_ratios':[0.435,0.13,0.435]})
    # 3 wavelengths, 5 angles
    
    modname = 'dirty_newforcebiasxi'
    xis = ['0.0','0.1','0.25','0.9','1.0']
    angles = ['000','090','180']

    custcmap = cubehelix.cmap(reverse = False, rot=1, start=0, 
                              startHue=320, sat = 2)

    wave = args.wave
    if wave == '035.11':
        minmax_y_vals = np.array([[0.1,1.0],
                                  [1e-2,1e0],
                                  [0.1,1.0]])
    elif wave == '151.99':
        minmax_y_vals = np.array([[100.,250.],
                                  [1e0,1e4],
                                  [100.,250.]])
    else:
        wave = '000.15'
        minmax_y_vals = np.array([[5e-4,1.5e-3],
                                  [1e-40,1.0],
                                  [1e-34,1e-31]])

    tau = '1e1'

    for i, xi in enumerate(xis):
        ifilenames = [modname + '/' + modname + '_' + xi +
                      '_slab_eff_t' + tau + '_i' + angle + 'a000_w' +
                      wave + '.fits' for angle in angles]

        for j, cfile in enumerate(ifilenames):
            # read in the image
            timage = np.squeeze(fits.getdata(cfile))

            if j == 1:
                timage = timage[90:240,90:510]
            else:
                timage = timage[90:510,90:510]

            #if timage.shape[0] == 300:
            #    timage = timage[9:255,45:257]
            #elif timage.shape[0] == 600:
            #    timage = timage[19:520,90:510]

            # display images
            tax = ax[j,i]
            cur_cax = tax.imshow(timage,
                                     norm=LogNorm(vmin=minmax_y_vals[j,0],
                                                  vmax=minmax_y_vals[j,1]),
                                     origin='lower',
                                     cmap=custcmap)
            tax.get_xaxis().set_visible(False)
            if j == 0:
                tax.set_title(r'$\xi_\mathrm{scat} = ' + xi + '$',
                                  fontsize=fontsize)
            if i == 0:
                tax.set_ylabel(r'$\theta = ' + angles[j] + '^\circ$')
                tax.yaxis.set_ticklabels([])
                tax.yaxis.set_ticks_position('none')
            else:
                tax.get_yaxis().set_visible(False)

    # add the overall label
    fig.text (0.05, 0.5, r'$\lambda = ' + wave + '$ $\mu m$', 
              horizontalalignment='center', rotation=90.,
              verticalalignment='center',fontsize=1.5*fontsize)
            
    # optimize the figure layout
    pyplot.tight_layout(rect=[0.06, 0, 1, 1])

    # display the plot
    save_name = modname+'_converge_t' + tau + '_w' + wave
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

