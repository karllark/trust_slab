#!/usr/bin/env python
#
# Code to make plots for comparing the global SEDs or images 
# for the TRUST BM1 Slab benchmark
#
# Written by: Karl Gordon (kgordon@stsci.edu)
#
# 23 Jul 2015 - written
# 26 Aug 2015 - merged the sed and image comparison codes into one
#               single code to avoid duplications
#
import os.path
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm

from astropy.io import fits

from comp_slab_models_sed_decompose import plot_decompose_sed
from comp_slab_models_image import plot_imagegrid

if __name__ == "__main__":

    good_angles = ['000','030','060','090','120','150','180','all']
    good_taus = ['1e-2','1e-1','1e0','1e1','all']
    good_waves = ['000.15','000.53','008.11','023.10','035.11','151.99','all']

    # commandline parser
    parser = argparse.ArgumentParser()

    # default is to compare the global SEDs, set for image comparison
    parser.add_argument("--image", help="Display image comparison", action="store_true")

    # options for which models to compare
    parser.add_argument("-w", "--wave", choices=good_waves, default='000.53',
                        help="wavelength to display")
    parser.add_argument("-t", "--tau", choices=good_taus, default='1e0', 
                        help="Optical depth of model run")
    parser.add_argument("-a", "--angle", choices=good_angles, default='090',
                        help="Viewing angle of model run")
    parser.add_argument("-m", "--max_plot_diff", metavar=float, default=10.0, 
                        help="maximum difference in percentage to plot")

    parser.add_argument("--eff", action="store_true",
                        help="Display results for the effective grain emission approximation [Default]")
    parser.add_argument("--equ", action="store_true",
                        help="Display results for the equilibrium heating only approximation")
    parser.add_argument("--sto", action="store_true",
                        help="Display results for the full heating solution (equilibrium + non-equilibrium)")

    parser.add_argument("--dirty_gtype", action="store_true",
                        help="display the different grain emission types (special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_econs", action="store_true",
                        help="energy conservation convergence (special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_mscat", action="store_true",
                        help="max scat convergence (special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_nphot", action="store_true",
                        help="nphot convergence (special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_forcebiasxi", action="store_true",
                        help="value of xi for biasing the forced scattering (special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_biasxi", action="store_true",
                        help="value of xi for biasing the regular scattering (special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_emitbiasxi", action="store_true",
                        help="value of xi for biasing the dust emission (special DIRTY runs) [default=False]")
    parser.add_argument("--dirty_nz", action="store_true",
                        help="number of z bins in slab (special DIRTY runs) [default=False]")
    parser.add_argument("--skirt_nz", action="store_true",
                        help="number of grid cells in the z direction (special SKIRT runs) [default=False]")
    parser.add_argument("--skirt_wr", action="store_true",
                        help="minimum weight reduction (special SKIRT runs) [default=False]")

    parser.add_argument("-e", "--eps", help="Save the plot as an encapsulated file",
                        action="store_true")
    parser.add_argument("-p", "--png", help="Save the plot as a portable network graphics file",
                        action="store_true")
    args = parser.parse_args()

    # check if 'all' then set to full set
    angles = [args.angle]
    if 'all' in angles:
        angles = good_angles[0:len(good_angles)-1]

    taus = [args.tau]
    if 'all' in taus:
        taus = good_taus[0:len(good_taus)-1]

    waves = [args.wave]
    if 'all' in waves:
        if args.sto:
            waves = ['000.15','000.53','008.11','023.10','151.99']
        else:
            waves = ['000.15','000.53','035.11','151.99']

    mplot_diff = float(args.max_plot_diff)

    # models to display
    plot_all = False 
    if args.dirty_nz:
        nbinzs = ['10','20','50','100','200','500']
        moddisplaynames = ['DI (Nz='+nbinz+')' for nbinz in reversed(nbinzs)]
        modnames = ['dirty_nbinz_'+nbinz for nbinz in reversed(nbinzs)]
        imodnames = ['dirty_nbinz/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'dirty_nz'
    elif args.skirt_nz:
        nbinzs = ['005','010','030','100','200','400']
        moddisplaynames = ['SK (Nz='+nbinz+')' for nbinz in reversed(nbinzs)]
        modnames = ['skirtnz'+nbinz for nbinz in reversed(nbinzs)]
        imodnames = ['skirt_nbinz/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'skirt_nz'
    elif args.dirty_nphot:
        nphots = ['1e5','3.2e5','1e6','3.2e6','1e7','3.2e7','1e8']
        moddisplaynames = ['DI (N='+nphot+')' for nphot in reversed(nphots)]
        modnames = ['dirty_nphot_'+nphot for nphot in reversed(nphots)]
        imodnames = ['dirty_nphot/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'dirty_nphot' 
    elif args.dirty_mscat:
        mscats = ['1','5','10','20','50','75','100','150','200','300','500','1000']
        moddisplaynames = ['DI (mscat=' + mscat + ')' for mscat in reversed(mscats)]
        modnames = ['dirty_mscat_' + mscat for mscat in reversed(mscats)]
        imodnames = ['dirty_mscat/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'dirty_mscat'
    elif args.skirt_wr:
        weightred = ['1e8','1e7','1e6','1e5','1e4','1e3']
        moddisplaynames = ['SK (wr='+wr+')' for wr in weightred]
        modnames = ['skirtwr'+wr for wr in weightred]
        imodnames = ['skirt_wr/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'skirt_wr'
    elif args.dirty_econs:
        econtargs = ['1.0','0.32','0.1','0.032','0.01','0.0032','0.001']
        moddisplaynames = ['DI (econs='+econtarg+')' for econtarg in reversed(econtargs)]
        modnames = ['dirty_econs_'+econtarg for econtarg in reversed(econtargs)]
        imodnames = ['dirty_econs/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'dirty_econs'
    elif args.dirty_forcebiasxi:
        xis = ['0.0','0.25','0.5','0.75','1.0']
        moddisplaynames = ['DI (fxis=' + xi + ')' for xi in xis]
        modnames = ['dirty_forcebiasxi_' + xi for xi in xis]
        imodnames = ['dirty_forcebiasxi/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'dirty_forcebiasxi'
    elif args.dirty_biasxi:
        xis = ['0.0','0.05','0.10','0.15','0.25']
        moddisplaynames = ['DI (xis=' + xi + ')' for xi in xis]
        modnames = ['dirty_biasxi_' + xi for xi in xis]
        imodnames = ['dirty_biasxi/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'dirty_forcebiasxi'
    elif args.dirty_emitbiasxi:
        xis = ['0.0','0.1','0.25','0.5','0.75','1.0']
        moddisplaynames = ['DI (exis=' + xi + ')' for xi in xis]
        modnames = ['dirty_emitbiasxi_' + xi for xi in xis]
        imodnames = ['dirty_emitbiasxi/' + modname + '_slab_eff' for modname in modnames]
        scomp = 0
        save_str = 'dirty_forcebiasxi'
    elif args.dirty_gtype:
        gtypes = ['equ','eff']
        moddisplaynames = ['DI (gtype='+gtype+')' for gtype in gtypes]
        modnames = ['dirty','dirty']
        imodnames = ['dirty/dirty_slab_' + gtype for gtype in gtypes]
        scomp = 0
        save_str = 'dirty_gtype'
        plot_all = True
    elif args.equ:
        moddisplaynames = ['CRT','DIRTY']
        modnames = ['crt','dirty']
        imodnames = [modname + '/' + modname + '_slab_equ' for modname in modnames]
        scomp = -1
        save_str = 'equ'
    elif args.sto:
        moddisplaynames = ['CRT']
        modnames = ['crt']
        imodnames = [modname + '/' + modname + '_slab_sto' for modname in modnames]
        scomp = -1
        save_str = 'sto'
    else:
        moddisplaynames = ['CRT','DART-ray','DIRTY','Hyperion','SKIRT','TRADING','SOC']
        modnames = ['crt','dartr','dirty','hyper','skirt','tradi','SOC']
        #moddisplaynames = ['CRT','DART-ray','DIRTY','DIRTY(old)','Hyperion','SKIRT','SOC','TRADING']
        #modnames = ['crt','dartr','dirty','dirty_prescat','hyper','skirt','SOC','tradi']
        #moddisplaynames = ['CRT','DART-ray','DIRTY','Hyperion','SKIRT','SOC','TRADING']
        #modnames = ['crt','dartr','dirty_prescat','hyper','skirt','SOC','tradi']
        imodnames = [modname + '/' + modname + '_slab_eff' for modname in modnames]
        #imodnames[3] = 'dirty_prescat/dirty_slab_eff'
        scomp = -1
        save_str = ''

    for angle in angles:
        for tau in taus:
            if args.image:
                for wave in waves:
                    plot_imagegrid(imodnames, moddisplaynames, wave, tau, angle,
                                   max_plot_diff=mplot_diff, comp_index=scomp,
                                   save_str=save_str, save_eps=args.eps, save_png=args.png)
            else:
                plot_decompose_sed(imodnames, moddisplaynames, tau, angle,
                                   single_comp=scomp, max_plot_diff=mplot_diff, plot_all=plot_all,
                                   save_str=save_str, save_eps=args.eps, save_png=args.png)
                                   
