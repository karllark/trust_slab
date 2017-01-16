#!/usr/bin/env python
#
# Code to plot the quantitiave comparisons of the global SEDs 
#          between different codes
#
# 27 Apr 2016: written (Karl Gordon)
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

from astropy.table import Table

def plot_indiv_comp(ax, tau, good_angles, tagstr, dispstr, compname, 
                    dcol, dsym='-', dustemis=False):

    n_angles = len(good_angles)
    vinit = False
    for k, angle in enumerate(good_angles):
        tab_name =  'slab_t' + tau + '_i' + angle + '_decomposed_sed_comp' + \
            tagstr
            
        cur_table = Table.read('dat/' + tab_name+'.dat',
                               format='ascii.commented_header')

        indxs, = np.where(cur_table['component'] == compname)

        if not vinit:
            n_models = len(indxs)
            modnames = np.array(cur_table['model'][indxs])
            mod_std = np.zeros((n_angles,n_models))
            mod_maxdev = np.zeros((n_angles,n_models))
            vinit = True
            
        mod_std[k,:] = cur_table['stddev'][indxs]
        mod_maxdev[k,:] = cur_table['maxabsdev'][indxs]

    not_used_models = ["DART-ray"]
    for k, model in enumerate(modnames):
        trans = 1.0
        if (model in not_used_models) and dustemis:
            trans = 0.25
        ax.plot(good_angles, mod_std[:,k], dcol[model]+dsym, 
                label=model+r' $\sigma$', alpha=trans)

    mvals = modnames
    indxs = np.arange(len(mvals))
    if dustemis:
        for nmod in not_used_models:
            nindxs, = np.where(mvals[indxs] != nmod)
            indxs = indxs[nindxs]

    print('Global SED: ', tau, compname, np.max(mod_std[:,indxs]))
            
    ax.set_title(r'$\tau_z =$ ' + tau + ' (' + dispstr + ', ' + compname + ')')
    ax.set_yscale('log')    
    ax.set_ylim(1e-3,1e4)
    ax.set_ylabel(r'$\bar{\Delta}$ [%]')
    ax.set_xlabel(r'$\theta$')

    return modnames

def plot_full_comp(compname, dcompnames, args, good_angles, ax,
                   tagstr, dispstr,
                   fontsize=14, dsym='-', plegend=True, dustemis=False):

    # setup the plot
    
    #taus = ['1e-2','1e-1','1e0','1e1']
    #fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(15,8))
    #tax = [ax[0,0],ax[0,1],ax[1,0],ax[1,1]]

    taus = ['1e0','1e1']
    tax = ax

    col = ['b','g','r','c','m','y','k']
    models = ['CRT','DART-ray','DIRTY','Hyperion','SKIRT',
              'TRADING','SOC']
    dcol = dict(zip(models, col))

    for k, cur_ax in enumerate(tax):
        pmodels = plot_indiv_comp(cur_ax, taus[k], good_angles, 
                                  tagstr, dispstr, 
                                  dcompnames[compname], dcol,
                                  dsym=dsym, dustemis=dustemis)
        if k == 0:
            save_pmodels = pmodels

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

    # make y-axis better
    tax[1].yaxis.tick_right()
    tax[1].yaxis.set_ticks_position('both')
    tax[1].yaxis.set_label_position("right")

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
                        help="Display results for the full heating solution" + \
                        " (equilibrium + non-equilibrium)")
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

    dcompnames = {}
    dcompnames['total'] = 'Total'
    dcompnames['dstel'] = 'Direct Stellar'
    dcompnames['dscat'] = 'Scattered Stellar'
    dcompnames['demis'] = 'Direct Dust Emission'
    dcompnames['demisscat'] = 'Scattered Dust Emission'
    dcompnames['tstel'] = 'Transparent Stellar'

    dispstr = 'eff'
    tagstr = ''
    if args.equ:
        dispstr = 'equ'
        tagstr = '_equ'
    elif args.sto:
        dispstr = 'sto'
        tagstr = '_sto'

    for cur_compname in dcompnames.keys():
        fig, ax = plt.subplots(ncols=2,figsize=(15,6))

        plot_full_comp(cur_compname, dcompnames, args, good_angles, 
                       ax, tagstr, dispstr)

        fig.tight_layout()

        save_name = 'slab_qcomp_' + dispstr + '_' + cur_compname
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

