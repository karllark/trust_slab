#!/usr/bin/env python
#
# Code to plot the quantitiave comparisons of the Y image slices
#          between different codes
#
# 17 Jun 2016: written (Karl Gordon)
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

from astropy.table import Table

def plot_indiv_slice(ax, tau, waves, angles, run_tag,
                     dcol, fontsize=16, irwaves=False):

    sym = ['-','--','-.',':']
    for n, wave in enumerate(waves):
        for k, angle in enumerate(angles):

            tab_name =  'dat/slab_t' + tau + '_i' + angle + '_w' + wave + \
                '_image_comp' + run_tag
            
            cur_table = Table.read(tab_name+'.dat',
                                   format='ascii.commented_header')

            # initialize and get model names
            if n == 0 and k == 0:
                mvals = np.array(cur_table['model'])
                pvals = np.empty((len(mvals),len(angles)))
            
            pvals[:,k] = np.array(cur_table['cut1_stddev'])

        not_used_models = ["DART-ray","TRADING"]
        for i, cmodel in enumerate(mvals):
            trans = 1.0
            if (tau == '1e1') and irwaves:
                if cmodel in not_used_models:
                    trans = 0.25
            ax.plot(angles, pvals[i,:], dcol[cmodel]+sym[n], alpha=trans,
                    label=cmodel + r' $\lambda = ' + wave + '$')

        indxs = np.arange(len(mvals))
        if (tau == '1e1') and irwaves:
            for nmod in not_used_models:
                nindxs, = np.where(mvals[indxs] != nmod)
                indxs = indxs[nindxs]

        print('Y slice: ', tau, wave, np.max(pvals[indxs,:]))
            
    #sones = np.full((len(angles)),1.0)
    #ax.plot(angles,sones,'k--')
    #ax.plot(angles,5.0*sones,'k-.')

    ax.set_yscale('log')    
    ax.set_ylim(0.5e-1,1e4)
    ax.set_ylabel(r'$\bar{\Delta}$ [%]')
    ax.set_xlabel(r'$\theta$')

    ax.set_title(r'$\tau_z =$ ' + tau + ' (Y slice)')

    return mvals

    #ax.legend(loc=3)

def plot_slice_all_taus(args, good_angles, waves, ax, dispstr,
                        fontsize=14, plegend=True, irwaves=False):

    # setup the plot
    taus = ['1e0','1e1']
    tax = ax

    col = ['b','g','r','c','m','y','k']
    models = ['CRT','DART-ray','DIRTY','Hyperion','SKIRT',
              'TRADING','SOC']
    dcol = dict(zip(models, col))

    for k, cur_ax in enumerate(tax):
        pmodels = plot_indiv_slice(cur_ax, taus[k], waves, good_angles, 
                                   dispstr,dcol, irwaves=irwaves)
        if k == 0:
            save_pmodels = pmodels

    # Create two custom legends (more compact)

    # models
    if plegend:
        arts = [plt.Line2D((0,1),(0,0), color=dcol[cmodel], linestyle='-') 
                for cmodel in save_pmodels]
        leg1 = tax[0].legend(arts,
                             save_pmodels,
                             fontsize=1.25*fontsize,
                             loc='upper left',
                             ncol=2)
        leg1.get_frame().set_linewidth(2)

    # Add the legend manually to the current Axes.
    #tax[0].add_artist(leg1)

    # waves
    if args.sto and not args.uvopt:
        arts = [plt.Line2D((0,1),(0,0), color='k', linestyle='-'),
                plt.Line2D((0,1),(0,0), color='k', linestyle='--'),
                plt.Line2D((0,1),(0,0), color='k', linestyle='-.')]
        labs = [r'$\lambda = '+waves[0]+'$',
                r'$\lambda = '+waves[1]+'$',
                r'$\lambda = '+waves[2]+'$']
    else:
        arts = [plt.Line2D((0,1),(0,0), color='k', linestyle='-'),
                plt.Line2D((0,1),(0,0), color='k', linestyle='--')]
        labs = [r'$\lambda = '+waves[0]+'$',
                r'$\lambda = '+waves[1]+'$']
        
    leg2 = tax[1].legend(arts, labs,
                         title='Y Slice',
                         fontsize=1.25*fontsize,
                         loc='upper center', 
                         bbox_to_anchor=(-0.05,0.97))
    leg2.get_title().set_fontsize(fontsize*1.5)
    leg2.get_frame().set_linewidth(2)

    # make y axis nice

    tax[1].yaxis.tick_right()
    tax[1].yaxis.set_ticks_position('both')
    tax[1].yaxis.set_label_position("right")

if __name__ == "__main__":

    good_angles = ['000','030','060','090','120','150','180']

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--uvopt", action="store_true",
                        help="Display UV/optical results")
    parser.add_argument("--eff", action="store_true",
                        help="Display results for the effective grain " + \
                        "emission approximation [Default]")
    parser.add_argument("--equ", action="store_true",
                        help="Display results for the equilibrium heating " + \
                        "only approximation")
    parser.add_argument("--sto", action="store_true",
                        help="Display results for the full heating solution " +  \
                        "(equilibrium + non-equilibrium)")
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

    # setup the plots
    fontsize = 14
    font = {'size'   : fontsize}

    mpl.rc('font', **font)

    mpl.rc('lines', linewidth=2)
    mpl.rc('axes', linewidth=2)
    mpl.rc('xtick.major', width=2)
    mpl.rc('ytick.major', width=2)

    if args.uvopt:
        waves = ['000.15','000.53']
        irwaves = False
    else:
        if args.sto:
            waves = ['008.11','023.10','151.99']
        else:
            waves = ['035.11','151.99']
        irwaves = True

    dispstr = 'eff'
    tagstr = ''
    if args.equ:
        dispstr = 'equ'
        tagstr = '_equ'
    elif args.sto:
        dispstr = 'sto'
        tagstr = '_sto'

    fig, ax = plt.subplots(ncols=2,figsize=(15,6))

    plot_slice_all_taus(args, good_angles, waves, ax, tagstr)

    fig.tight_layout()

    save_name = 'slab_qcomp_image_' + dispstr

    if args.uvopt:
        save_name += '_uvopt'
    else:
        save_name += '_ir'

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
