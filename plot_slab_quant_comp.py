#!/usr/bin/env python
#
# Code to plot the quantitiave comparisons between different codes
#
# 27 Apr 2016: written (Karl Gordon)
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import matplotlib as mpl

from astropy.table import Table

def plot_indiv_comp(ax, tau, good_angles, tagstr, dispstr, compname):

    n_angles = len(good_angles)
    vinit = False
    for k, angle in enumerate(good_angles):
        tab_name =  'slab_t' + tau + '_i' + angle + '_decomposed_sed_comp' + tagstr
            
        cur_table = Table.read(tab_name+'.dat',
                               format='ascii.commented_header')

        indxs, = np.where(cur_table['component'] == compname)

        if not vinit:
            n_models = len(indxs)
            modnames = np.array(cur_table['model'][indxs])
            mod_std = np.zeros((n_angles,n_models))
            vinit = True
            
        mod_std[k,:] = cur_table['stddev'][indxs]

    for k, model in enumerate(modnames):

        ax.plot(good_angles, mod_std[:,k], label=model)
            

    #ax.plot([min_val,max_val],[1.0,1.0],'k--')
    #ax.plot([min_val,max_val],[5.0,5.0],'k-.')

    ax.set_title(r'$\tau =$ ' + tau + ' (' + dispstr + ', ' + compname + ')')
    ax.set_yscale('log')    
    ax.set_ylim(1e-3,1e2)
    ax.set_ylabel(r'$\sigma$ [%]')
    ax.set_xlabel('angle')


def plot_full_comp(compname, dcompnames, args, good_angles):

    # setup the plot
    fig, ax = pyplot.subplots(nrows=2,ncols=2,figsize=(15,10))
    
    taus = ['1e-2','1e-1','1e0','1e1']
    tax = [ax[0,0],ax[0,1],ax[1,0],ax[1,1]]

    dispstr = 'eff'
    tagstr = ''
    if args.equ:
        dispstr = 'equ'
        tagstr = '_equ'
    elif args.sto:
        dispstr = 'sto'
        tagstr = '_sto'

    for k, cur_ax in enumerate(tax):
        plot_indiv_comp(cur_ax, taus[k], good_angles, tagstr, dispstr, 
                        dcompnames[compname])

    tax[0].legend(loc=2, ncol=2)

    fig.tight_layout()


    save_name = 'slab_qcomp_' + dispstr + '_' + compname
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

if __name__ == "__main__":

    good_angles = ['000','030','060','090','120','150','180']

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--eff", action="store_true",
                        help="Display results for the effective grain " + \
                        "emission approximation [Default]")
    parser.add_argument("--equ", action="store_true",
                        help="Display results for the equilibrium heating " + \
                        "only approximation")
    parser.add_argument("--sto", action="store_true",
                        help="Display results for the full heating solution " + \
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
    fontsize = 12
    font = {'size'   : fontsize}

    mpl.rc('font', **font)

    mpl.rc('lines', linewidth=2)
    mpl.rc('axes', linewidth=2)
    mpl.rc('xtick.major', width=2)
    mpl.rc('ytick.major', width=2)

    dcompnames = {}
    dcompnames['total'] = 'Total'
    dcompnames['dstel'] = 'Direct Stellar'
    dcompnames['dscat'] = 'Scattered Stellar'
    dcompnames['demis'] = 'Direct Dust Emission'
    dcompnames['demisscat'] = 'Scattered Dust Emiss'
    dcompnames['tstel'] = 'Transparent Stellar'

    for cur_compname in dcompnames.keys():
        plot_full_comp(cur_compname, dcompnames, args, good_angles)
