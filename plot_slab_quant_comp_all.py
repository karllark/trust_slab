#!/usr/bin/env python
#
# Code to plot the quantitiave comparisons of the global SEDs and
#          Y image slices between different codes
#
# 17 Jun 2016: written (Karl Gordon)
# 25 Oct 2017: updated to plot all 6 panels
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

from astropy.table import Table

from plot_slab_quant_comp_image import plot_slice_all_taus
from plot_slab_quant_comp import plot_full_comp

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

    dispstr = 'eff'
    tagstr = ''
    if args.equ:
        dispstr = 'equ'
        tagstr = '_equ'
    elif args.sto:
        dispstr = 'sto'
        tagstr = '_sto'

    fig, ax = plt.subplots(ncols=2, nrows=2,figsize=(15,12))

    # plot the global SED comparisons

    dcompnames = {}
    dcompnames['total'] = 'Total'
    dcompnames['dstel'] = 'Direct Stellar'
    dcompnames['dscat'] = 'Scattered Stellar'
    dcompnames['demis'] = 'Direct Dust Emission'
    dcompnames['demisscat'] = 'Scattered Dust Emiss'
    dcompnames['tstel'] = 'Transparent Stellar'

    if args.uvopt:
        compnames = ['dscat','dstel']
        dustemis = False
    else:
        compnames = ['demis','demisscat']
        dustemis = True
        
    syms = ['-','--']
    leg = [False, False]
    for k, cur_compname in enumerate(compnames):
        plot_full_comp(cur_compname, dcompnames, args, good_angles, 
                       [ax[0,0],ax[0,1]], tagstr, dispstr,
                       fontsize=fontsize, dsym=syms[k], plegend=leg[k],
                       dustemis=dustemis)

    # remove x axes and make set plot titles
    tau = ['1','10']
    for i in range(2):
        ax[0,i].set_xlabel('')
        ax[0,i].xaxis.set_ticklabels([])
        ax[0,i].set_title(r'$\tau_z =$ ' + tau[i])

    # Add the legend manually to distingish between the two global SED results
    dcompnames['demisscat'] = 'Scattered Dust Emission'
    arts = [plt.Line2D((0,1),(0,0), color='k', linestyle=csym) for csym in syms]
    labs = [dcompnames[cur_compname] for cur_compname in compnames]

    if args.uvopt:
        ltitle = 'Global SED'
    else:
        ltitle = 'Global SED (' + dispstr + ' case)'

    tax = ax[0,1]
    leg2 = tax.legend(arts, labs,
                      title=ltitle,
                      fontsize=1.25*fontsize,
                      loc='upper center', 
                      bbox_to_anchor=(-0.05,0.97))
    leg2.get_title().set_fontsize(fontsize*1.5)
    leg2.get_frame().set_linewidth(2)

    # plot the Y slice comparisons

    if args.uvopt:
        waves = ['000.15','000.53']
        irwaves = False
    else:
        if args.sto:
            waves = ['008.11','023.10','151.99']
        else:
            waves = ['035.11','151.99']
        irwaves = True

    plot_slice_all_taus(args, good_angles, waves, 
                        [ax[1,0],ax[1,1]], tagstr,
                        fontsize=fontsize, 
                        plegend=True, irwaves=irwaves)

    # remove plot titles
    for i in range(2):
        ax[1,i].set_title('')

    fig.tight_layout()

    save_name = 'slab_qcomp_image_all_' + dispstr

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
