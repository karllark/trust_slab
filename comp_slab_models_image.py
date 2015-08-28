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
                   comp_index=-2, max_plot_diff=100.0, save_str='',
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
    symtype = ['b-','g-','r-','c-','m-','y-','k-','b--','g--','r--','c--','m--','y--','k--']
    fontsize = 12

    # cut info
    cut1 = np.array([90,110])
    cut2 = np.array([80,100])
    
    # decide the number of columns for the images
    nrows = 4
    n_image_col = 2
    dm = divmod(n_files,nrows)
    if dm[0] >= n_image_col & dm[0] > 0:
        n_image_col = dm[0] + 1

    # setup figure
    xsize = 12
    ysize = 12
    if n_image_col > 2:
        xsize += 3.0*(n_image_col-2)
    fig, ax = pyplot.subplots(figsize=(xsize,ysize))

    # use gridspec to allow for one plot to be larger than the others
    gs = gridspec.GridSpec(nrows, n_image_col+3, width_ratios=2*[1.5] + n_image_col*[1.0] + [0.15])
    ax = []
    for i in range(n_files):
        gs_indxs = np.array(divmod(fileindxs[i],n_image_col))
        gs_indxs[1] += 2
        ax.append(pyplot.subplot(gs[gs_indxs[0],gs_indxs[1]]))
    # cut plots
    ax.append(pyplot.subplot(gs[0,0:2]))
    ax.append(pyplot.subplot(gs[1,0:2]))
    ax.append(pyplot.subplot(gs[2,0:2]))
    ax.append(pyplot.subplot(gs[3,0:2]))
        
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
            nx = timage.shape[0]
            ny = timage.shape[1]
            all_images = np.empty((timage.shape[0],timage.shape[1],n_files))
            minmax_vals = np.empty((n_files,2))
            if timage.shape[0] > 300:
                cut1 *= 2
                cut2 *= 2

        all_images[:,:,i] = timage

    # get the image for comparison (if desired)
    ave_image_comp = np.median(all_images,axis=2)
    if comp_index > -2:
        if comp_index >= 0:
            ave_image_comp = np.array(all_images[:,:,comp_index])
        #for i in range(n_files):
        #    all_images[:,:,i] = all_images[:,:,i] - ave_image_comp

        #    timage = all_images[:,:,i]
        #    gindxs = np.where(ave_image_comp[:] > 0.)
        #    timage[gindxs] /= ave_image_comp[gindxs]

    # get the min/max to plot
    for i in range(n_files):
        timage = all_images[:,:,i]
        #if comp_index > -2:
        #    gindxs = np.where(np.isfinite(timage[:]))
        #else:
        gindxs = np.where((timage[:] > 0.) & (timage[:] < np.max(timage[:])))
        if len(gindxs[0]) > 0:
            minmax_vals[i,0] = np.min(timage[gindxs])
            minmax_vals[i,1] = np.max(timage[gindxs])
        else:
            print(i,cfile,' has no positive values')

    #if comp_index > -2:
    #    plot_minmax = [-0.1,0.1]
    #else:
    gindxs, = np.where(minmax_vals[:,0] > 0)
    plot_minmax = [np.median(minmax_vals[gindxs,0]),np.median(minmax_vals[gindxs,1])]

    # cut plot variables
    cut1_minmax_x_vals = np.array([1e6,-1e6])
    cut1_minmax_y_vals = np.array([1e6,-1e6])
    cut1_plot_x = np.arange(nx)

    cut2_minmax_x_vals = np.array([1e6,-1e6])
    cut2_minmax_y_vals = np.array([1e6,-1e6])
    cut2_plot_x = np.arange(ny)

    cut1_plot_y_all = np.median(ave_image_comp[:,cut1[0]:cut1[1]],axis=1)
    cut2_plot_y_all = np.median(ave_image_comp[cut2[0]:cut2[1],:],axis=0)

    # now display the cuts    
    for i in range(n_files):

        # first cut (y)
        cut1_plot_y = np.median(all_images[:,cut1[0]:cut1[1],i],axis=1)
        gindxs, = np.where(cut1_plot_y > 0.)
        if gindxs.shape[0] > 0:
            gindxs2 = gindxs[1:-1]
            cut1_minmax_x_vals[0] = np.min([cut1_minmax_x_vals[0],np.min(cut1_plot_x[gindxs2])])
            cut1_minmax_x_vals[1] = np.max([cut1_minmax_x_vals[1],np.max(cut1_plot_x[gindxs2])])
            cut1_minmax_y_vals[0] = np.min([cut1_minmax_y_vals[0],np.min(cut1_plot_y[gindxs2])])
            cut1_minmax_y_vals[1] = np.max([cut1_minmax_y_vals[1],np.max(cut1_plot_y[gindxs2])])
            if i < 6:
                tname = displaynames[i]
            else:
                tname = None
            ax[n_files].plot(cut1_plot_x[gindxs],cut1_plot_y[gindxs],symtype[fileindxs[i]],label=tname)

            # percent difference plot for first cut
            gindxs, = np.where((cut1_plot_y > 0) & (cut1_plot_y_all > 0))
            if i >= 6:
                tname = displaynames[i]
            else:
                tname = None
            ax[n_files+1].plot(cut1_plot_x[gindxs],100.*(cut1_plot_y[gindxs]-cut1_plot_y_all[gindxs])/cut1_plot_y_all[gindxs],symtype[fileindxs[i]],label=tname)
            
        # second cut (x)
        cut2_plot_y = np.median(all_images[cut2[0]:cut2[1],:,i],axis=0)
        gindxs, = np.where(cut2_plot_y > 0.)
        if gindxs.shape[0] > 0:
            gindxs2 = gindxs[1:-1]
            cut2_minmax_x_vals[0] = np.min([cut2_minmax_x_vals[0],np.min(cut2_plot_x[gindxs2])])
            cut2_minmax_x_vals[1] = np.max([cut2_minmax_x_vals[1],np.max(cut2_plot_x[gindxs2])])
            cut2_minmax_y_vals[0] = np.min([cut2_minmax_y_vals[0],np.min(cut2_plot_y[gindxs2])])
            cut2_minmax_y_vals[1] = np.max([cut2_minmax_y_vals[1],np.max(cut2_plot_y[gindxs2])])
            if i < 6:
                tname = displaynames[i]
            else:
                tname = None
            ax[n_files+2].plot(cut2_plot_x[gindxs],cut2_plot_y[gindxs],symtype[fileindxs[i]],label=tname)
        
            # percent difference plot for first cut
            gindxs, = np.where((cut2_plot_y > 0) & (cut2_plot_y_all > 0))
            if i >= 6:
                tname = displaynames[i]
            else:
                tname = None
            ax[n_files+3].plot(cut2_plot_x[gindxs],100.*(cut2_plot_y[gindxs]-cut2_plot_y_all[gindxs])/cut2_plot_y_all[gindxs],symtype[fileindxs[i]],label=tname)

    # setup for the first cut plot
    ax[n_files].set_yscale('log')
    cut1_minmax_x_vals[0] -= 0.1*(cut1_minmax_x_vals[1] - cut1_minmax_x_vals[0])
    cut1_minmax_x_vals[1] += 0.5*(cut1_minmax_x_vals[1] - cut1_minmax_x_vals[0])
    ax[n_files].set_xlim(cut1_minmax_x_vals)
    cut1_minmax_y_vals[0] = 10**(np.log10(cut1_minmax_y_vals[0]) - (0.1*(np.log10(cut1_minmax_y_vals[1]) - np.log10(cut1_minmax_y_vals[0]))))
    cut1_minmax_y_vals[1] = 10**(np.log10(cut1_minmax_y_vals[1]) + (0.1*(np.log10(cut1_minmax_y_vals[1]) - np.log10(cut1_minmax_y_vals[0]))))
    ax[n_files].set_ylim(cut1_minmax_y_vals)
    ax[n_files].set_ylabel('SB [MJy/sr]')
    ax[n_files].set_title('Y slice ($'+str(cut1[0])+' \leq x \leq '+str(cut1[1])+ '$)')
    ax[n_files].legend(loc=1,fontsize=fontsize)

    ax[n_files+1].set_ylabel('% difference')
    ax[n_files+1].set_xlim(cut1_minmax_x_vals)
    cur_ylim = ax[n_files+1].get_ylim()
    new_ylim = [max([cur_ylim[0],-1.0*max_plot_diff]),min([cur_ylim[1],max_plot_diff])]
    ax[n_files+1].set_ylim(new_ylim)
    if n_files > 6:
        ax[n_files+1].legend(loc=1,fontsize=fontsize)

    # setup for the second cut plot
    ax[n_files+2].set_yscale('log')
    cut2_minmax_x_vals[0] -= 0.1*(cut2_minmax_x_vals[1] - cut2_minmax_x_vals[0])
    cut2_minmax_x_vals[1] += 0.5*(cut2_minmax_x_vals[1] - cut2_minmax_x_vals[0])
    ax[n_files+2].set_xlim(cut2_minmax_x_vals)
    cut2_minmax_y_vals[0] = 10**(np.log10(cut2_minmax_y_vals[0]) - (0.1*(np.log10(cut2_minmax_y_vals[1]) - np.log10(cut2_minmax_y_vals[0]))))
    cut2_minmax_y_vals[1] = 10**(np.log10(cut2_minmax_y_vals[1]) + (0.1*(np.log10(cut2_minmax_y_vals[1]) - np.log10(cut2_minmax_y_vals[0]))))
    ax[n_files+2].set_ylim(cut2_minmax_y_vals)
    ax[n_files+2].set_ylabel('SB [MJy/sr]')
    ax[n_files+2].set_title('X slice ($'+str(cut2[0])+' \leq y \leq '+str(cut2[1])+ '$)')
    ax[n_files+2].legend(loc=1,fontsize=fontsize)

    ax[n_files+3].set_ylabel('% difference')
    ax[n_files+3].set_xlim(cut2_minmax_x_vals)
    cur_ylim = ax[n_files+1].get_ylim()
    new_ylim = [max([cur_ylim[0],-1.0*max_plot_diff]),min([cur_ylim[1],max_plot_diff])]
    ax[n_files+3].set_ylim(new_ylim)
    if n_files > 6:
        ax[n_files+3].legend(loc=1,fontsize=fontsize)

    # now display the images    
    for i in range(n_files):
        # display images
        #if comp_index > -2:
        #    cur_cax = ax[i].imshow(all_images[:,:,i],vmin=plot_minmax[0],vmax=plot_minmax[1], origin='lower')#,
        #else:
        #cur_cax = ax[i].imshow(all_images[:,:,i],norm=LogNorm(vmin=plot_minmax[0],vmax=plot_minmax[1]), origin='lower')#,
        cur_cax = ax[i].imshow(all_images[:,:,i],norm=LogNorm(vmin=cut1_minmax_y_vals[0],vmax=cut1_minmax_y_vals[1]), origin='lower')#,
#                               cmap=pyplot.get_cmap('cubehelix'))
        ax[i].set_title(displaynames[i],fontsize=fontsize)
        ax[i].get_xaxis().set_visible(False)
        ax[i].get_yaxis().set_visible(False)


    # add the overall label
    fig.text (0.5, 0.99, fig_label, horizontalalignment='center',
              verticalalignment='top',fontsize=1.5*fontsize)

    # colorbar
    fig.colorbar(cur_cax, cax=(pyplot.subplot(gs[0:nrows,n_image_col+2])))
    
    # optimize the figure layout
    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.96])

    # display the plot
    save_name =  'slab_t' + tau + '_i' + angle + '_w' + wave + '_image_comp'

    if save_str != '':
        save_name += '_' + save_str
    
    if save_png:
        fig.savefig(save_name+'.png')
        fig.savefig(save_name+'_small.png',dpi=13)
    elif save_eps:
        fig.savefig(save_name+'.eps')
    else:
        pyplot.show()

    pyplot.close(fig)

if __name__ == "__main__":

    print('Use comp_slab_models.py --image')
