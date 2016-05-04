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

def plot_indiv_comp(ax, taus, angles, compname):

    for tau in taus:
        for angle in angles:
            tab_name =  'slab_t' + tau + '_i' + angle + \
                '_decomposed_sed_comp_scomp0_' + run_tag
            
            cur_table = Table.read(tab_name+'.dat',
                                   format='ascii.commented_header')

            indxs, = np.where(cur_table['component'] == compname)

        # extract the number of bins from the model name
            mvals = np.empty((len(indxs)))
            for i, cmodel in enumerate(cur_table['model'][indxs]):
                eq_pos = cmodel.find('=')
                pa_pos = cmodel.find(')')
                mvals[i] = float(cmodel[eq_pos+1:pa_pos])

            nvals = len(mvals)
            ax.plot(mvals[1:nvals], cur_table['stddev'][indxs[1:nvals]], 
                    label=r'$\tau = ' + tau + '; i = ' + angle + '$')
            
    min_val = min(mvals[1:nvals])
    max_val = max(mvals[1:nvals])

    ax.plot([min_val,max_val],[1.0,1.0],'k--')
    ax.plot([min_val,max_val],[5.0,5.0],'k-.')

    ax.set_title(compname)
    ax.set_xscale('log')    
    ax.set_yscale('log')    
    #ax.set_ylim(0.5e-1,1e2)
    ax.set_ylabel(r'$\sigma$ [%]')
    ax.set_xlabel(kxlabel)
    ax.legend()


if __name__ == "__main__":

    good_angles = ['000','090','180']

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--dirty_nz", action="store_true",
                        help="number of z bins in slab " + \
                        "(special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_nxy", action="store_true",
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
    fig, ax = pyplot.subplots(nrows=2,ncols=2,figsize=(15,10))
    
    # setup the plots
    fontsize = 12
    font = {'size'   : fontsize}

    mpl.rc('font', **font)

    mpl.rc('lines', linewidth=2)
    mpl.rc('axes', linewidth=2)
    mpl.rc('xtick.major', width=2)
    mpl.rc('ytick.major', width=2)

    # read in the table data for each
    if args.dirty_nz:
        taus = ['1e0','1e1']
        angles = ['090'] 
        run_tag = 'dirty_nz'
        kxlabel = r'$n_z$'
    elif args.dirty_nxy:
        taus = ['1e0','1e1']
        angles = ['000','180']
        #angles = ['000']
        run_tag = 'dirty_nxy'
        kxlabel = r'$n_{xy}$'
    else:  # do the # photons case
        taus = ['1e0','1e1']
        #taus = ['1e0']
        angles = ['000','090','180']
        run_tag = 'dirty_nphot'
        kxlabel = r'$n_p$'

    dcompnames = {}
    dcompnames['total'] = 'Total'
    dcompnames['dscat'] = 'Scattered Stellar'
    dcompnames['demis'] = 'Direct Dust Emission'
    dcompnames['demisscat'] = 'Scattered Dust Emiss'

    tax = [ax[0,0],ax[0,1],ax[1,0],ax[1,1]]
    for k, cur_compname in enumerate(['total','dscat','demis','demisscat']):

        plot_indiv_comp(tax[k], taus, angles, dcompnames[cur_compname])

    fig.tight_layout()

    save_name = 'slab_converge_sed_' + run_tag
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