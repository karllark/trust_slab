#!/usr/bin/env python
#
# Code to plot the results of the convergence tests
#
# 26 Feb 2016: written (Karl Gordon)
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import matplotlib as mpl

from astropy.table import Table

def plot_converge_slice(ax, taus, waves, angle,
                        run_tag, kxlabel, plot_xlog=True,
                        fontsize=16):

    col = ['c','m','y','k']
    lstyle = ['-','--']
    for m, tau in enumerate(taus):
        for n, wave in enumerate(waves):
            tab_name =  'dat/slab_t' + tau + '_i' + angle + '_w' + wave + \
                        '_image_comp_' + run_tag
            
            cur_table = Table.read(tab_name+'.dat',
                                   format='ascii.commented_header')

            # extract the number of bins from the model name
            mvals = np.empty((len(cur_table['model'])))
            for i, cmodel in enumerate(cur_table['model']):
                eq_pos = cmodel.find('=')
                pa_pos = cmodel.find(')')
                mvals[i] = float(cmodel[eq_pos+1:pa_pos])

            nvals = len(mvals)
            ax.plot(mvals[1:nvals], cur_table['cut1_stddev'][1:nvals], 
                    col[n]+lstyle[m],
                    label=r'$\tau_z = ' + tau + '$; $\lambda = ' + wave + '$')
            
    min_val = min(mvals[1:nvals])
    max_val = max(mvals[1:nvals])

    ax.plot([min_val,max_val],[1.0,1.0],'k--')
    ax.plot([min_val,max_val],[5.0,5.0],'k-.')

    if plot_xlog:
        ax.set_xscale('log')    
    ax.set_yscale('log')    
    ax.set_ylim(0.5e-1,1e2)
    ax.set_ylabel(r'$\sigma$ [%]')
    ax.set_xlabel(kxlabel)
    #ax.legend(loc=3)


if __name__ == "__main__":

    good_angles = ['000','090','180']

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--angle", choices=good_angles, default='090',
                        help="Viewing angle of model run")
    parser.add_argument("--dirty_nz", action="store_true",
                        help="number of z bins in slab " + \
                        "(special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_nxy", action="store_true",
                        help="number of xy bins in slab " + \
                        "(special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_biasxi", action="store_true",
                        help="number of xy bins in slab " + \
                        "(special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_emitbiasxi", action="store_true",
                        help="number of xy bins in slab " + \
                        "(special DIRTY runs) [default=False]")
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

    # setup the plot
    fig, ax = pyplot.subplots(figsize=(12,10))
    
    # setup the plots
    fontsize = 16
    font = {'size'   : fontsize}

    mpl.rc('font', **font)

    mpl.rc('lines', linewidth=2)
    mpl.rc('axes', linewidth=2)
    mpl.rc('xtick.major', width=2)
    mpl.rc('ytick.major', width=2)

    # read in the table data for each
    angle = args.angle
    plot_xlog = True
    if args.dirty_nz:
        taus = ['1e0','1e1']
        waves = ['035.11','151.99']
        run_tag = 'dirty_nz'
        kxlabel = r'$n_z$'
    elif args.dirty_nxy:
        taus = ['1e0','1e1']
        waves = ['035.11','151.99']
        run_tag = 'dirty_nxy'
        kxlabel = r'$n_{xy}$'
    elif args.dirty_biasxi:
        taus = ['1e0','1e1']
        waves = ['000.15','000.53','035.11','151.99']
        run_tag = 'dirty_newforcebiasxi'
        kxlabel = r'$\xi_\mathrm{scat}$'
        plot_xlog = False
    elif args.dirty_emitbiasxi:
        taus = ['1e0','1e1']
        waves = ['000.15','000.53','035.11','151.99']
        run_tag = 'dirty_emitbiasxi'
        kxlabel = r'$\xi_\mathrm{emit}$'
        plot_xlog = False
    else:  # do the # photons case
        taus = ['1e0','1e1']
        #taus = ['1e0']
        waves = ['000.15','000.53','035.11','151.99']
        run_tag = 'dirty_nphot'
        kxlabel = r'$n_p$'

    plot_converge_slice(ax, taus, waves, angle, run_tag, kxlabel, plot_xlog=plot_xlog)

    if plot_xlog:
        ax.legend(loc=3)
    else:
        ax.legend(loc='upper center')
    fig.tight_layout()

    save_name = 'slab_converge_' + run_tag + '_i' + angle
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
