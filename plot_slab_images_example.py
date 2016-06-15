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
    fig, ax = pyplot.subplots(nrows=2,ncols=7,figsize=(15,5))
    # 3 wavelengths, 5 angles
    
    modname = 'skirt'
    waves = ['000.15','151.99']
    angles = ['000','030','060','090','120','150','180']

    custcmap = cubehelix.cmap(reverse = False, rot=1, start=0, 
                              startHue=320, sat = 2)

    minmax_y_vals = np.array([[1e-4,1e-1],
                              [1e0,3e2]])
    for i, wave in enumerate(waves):
        tau = '1e-1'
        ifilenames = [modname + '/' + modname + '_slab_sto_t' + tau + '_i' + 
                      angle + 'a000_w' +
                      wave + '.fits' for angle in angles]

        for j, cfile in enumerate(ifilenames):
            # read in the image
            timage = np.squeeze(fits.getdata(cfile))

            timage = timage[9:255,45:257]
            #if timage.shape[0] == 300:
            #elif timage.shape[0] == 600:
            #    timage = timage[19:520,90:510]

            # display images
            cur_cax = ax[i,j].imshow(timage,
                                     norm=LogNorm(vmin=minmax_y_vals[i,0],
                                                  vmax=minmax_y_vals[i,1]),
                                     origin='lower',
                                     cmap=custcmap)
            ax[i,j].get_xaxis().set_visible(False)
            if i == 0:
                ax[i,j].set_title(r'$\theta = ' + angles[j] + '^\circ$',
                                  fontsize=fontsize)
            if j == 0:
                ax[i,j].set_ylabel(r'$\lambda = ' + wave + '$ $\mu m$')
                ax[i,j].yaxis.set_ticklabels([])
                ax[i,j].yaxis.set_ticks_position('none')
            else:
                ax[i,j].get_yaxis().set_visible(False)

            
    # optimize the figure layout
    pyplot.tight_layout()

    # display the plot
    save_name = modname+'_images_example'
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

