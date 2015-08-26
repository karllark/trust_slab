#!/usr/bin/env python
#
# Code to make plots for comparing the images
# for the TRUST BM1 Slab benchmark
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 23 Jul 2015 - written
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm

from astropy.io import fits

def plot_imagegrid(modnames, moddisplaynames, wave, tau, angle,
                   comp_index=-1, 
                   save_eps=False, save_png=False):

    # generate the filename
    ifilenames = [modname + '_t' + tau + '_i' + angle + 'a000_w' + wave + '.fits'
                 for modname in modnames]
    n_orig_files = len(ifilenames)

    # check all the files exisit, adjust if not
    filenames = []
    displaynames = []
    fileindxs = []
    for i, cfile in enumerate(ifilenames):
        if os.path.isfile(cfile):
            filenames.append(cfile)
            displaynames.append(moddisplaynames[i])
            fileindxs.append(i)
    n_files = len(filenames)

    if n_files == 0:
        print('no files present')
        print(ifilenames)
        exit(0)
    
    # plot information
    fig_label = r'Slab, $\tau (1 \mu m)$ = '+tau+r', $\theta$ = ' + angle + ', $\lambda$ = ' + wave
    symtype = ['b-','g-','r-','c-','m-','y-','k-']
    total_symtype = ['k-','b--','b:','r--','r:','b--']
    fontsize = 12

    # setup figure
    fig, ax = pyplot.subplots(figsize=(10,10))

    # use gridspec to allow for one plot to be larger than the others
    # may do this later
    dm = divmod(n_files,3)
    if dm[0] >= 3 & dm[1] > 0:
        nrows = dm[0] + 1
    else:
        nrows = 3
    gs = gridspec.GridSpec(nrows, 4, width_ratios=[1.,1.,1.0,0.15])
    ax = []
    for i in range(n_files):
        ax.append(pyplot.subplot(gs[divmod(fileindxs[i],3)]))

    # read in the results from each model
    for i, cfile in enumerate(filenames):

        # read in the image
        timage = np.squeeze(fits.getdata(cfile))

        if displaynames[i] in ['CRT','SOC']:
            timage = np.rot90(timage,3)

        if timage.shape[0] == 300:
            timage = timage[9:260,45:255]
        elif timage.shape[0] == 600:
            timage = timage[19:520,90:510]

        # save the image in a cube
        if i == 0:
            all_images = np.empty((timage.shape[0],timage.shape[1],n_files))
            minmax_vals = np.empty((n_files,2))
        all_images[:,:,i] = timage

    # get the image for comparison (if desired)
    if comp_index > -2:
        if comp_index >= 0:
            ave_image_comp = np.array(all_images[:,:,comp_index])
        else:
            ave_image_comp = np.median(all_images,axis=2)
        for i in range(n_files):
            all_images[:,:,i] = all_images[:,:,i] - ave_image_comp

            timage = all_images[:,:,i]
            gindxs = np.where(ave_image_comp[:] > 0.)
            timage[gindxs] /= ave_image_comp[gindxs]

    # get the min/max to plot
    for i in range(n_files):
        timage = all_images[:,:,i]
        if comp_index > -2:
            gindxs = np.where(np.isfinite(timage[:]))
        else:
            gindxs = np.where((timage[:] > 0.) & (timage[:] < np.max(timage[:])))
        if len(gindxs[0]) > 0:
            minmax_vals[i,0] = np.min(timage[gindxs])
            minmax_vals[i,1] = np.max(timage[gindxs])
        else:
            print(i,cfile,' has no positive values')

    if comp_index > -2:
        plot_minmax = [-0.1,0.1]
    else:
        gindxs, = np.where(minmax_vals[:,0] > 0)
        plot_minmax = [np.median(minmax_vals[gindxs,0]),np.median(minmax_vals[gindxs,1])]
    
    for i in range(n_files):
        if comp_index > -2:
            cur_cax = ax[i].imshow(all_images[:,:,i],vmin=plot_minmax[0],vmax=plot_minmax[1], origin='lower')#,
        else:
            cur_cax = ax[i].imshow(all_images[:,:,i],norm=LogNorm(vmin=plot_minmax[0],vmax=plot_minmax[1]), origin='lower')#,
#                               cmap=pyplot.get_cmap('cubehelix'))
        ax[i].set_title(displaynames[i],fontsize=fontsize)
        ax[i].get_xaxis().set_visible(False)
        ax[i].get_yaxis().set_visible(False)

    # add the overall label
    fig.text (0.5, 0.99, fig_label, horizontalalignment='center',
              verticalalignment='top',fontsize=1.5*fontsize)

    # colorbar
    fig.colorbar(cur_cax, cax=(pyplot.subplot(gs[0:nrows,3])))
    
    # optimize the figure layout
    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.96])

    # display the plot
    save_name =  'slab_t' + tau + '_i' + angle + '_w' + wave + '_image_comp'
    
    if save_png:
        fig.savefig(save_name+'.png')
        fig.savefig(save_name+'_small.png',dpi=13)
    elif save_eps:
        fig.savefig(save_name+'.eps')
    else:
        pyplot.show()

    pyplot.close(fig)

if __name__ == "__main__":

    print('Use comp_slab_models.py')
