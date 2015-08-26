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
    good_waves = ['000.15','000.53','035.11','151.99','all']

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

    parser.add_argument("--econs", action="store_true",
                        help="energy conservation convergence (special DIRTY runs) [default=False]")
    parser.add_argument("--mscat", action="store_true",
                        help="max scat convergence (special DIRTY runs) [default=False]")
    parser.add_argument("--nphot", action="store_true",
                        help="nphot convergence (special DIRTY runs) [default=False]")
    parser.add_argument("--stau", action="store_true",
                        help="subdivision tau for clumps (special DIRTY runs) [default=False]")
    parser.add_argument("--nbinz", action="store_true",
                        help="number of z bins in slab (special DIRTY runs) [default=False]")
    parser.add_argument("--nz", action="store_true",
                        help="number of grid cells in the z direction (special SKIRT runs) [default=False]")
    parser.add_argument("--wr", action="store_true",
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
        waves = good_waves[0:len(good_waves)-1]

    mplot_diff = float(args.max_plot_diff)

    # models to display
    if args.stau:
        moddisplaynames = ['DIRTY (Nz=400)','DIRTY (Nz=200)','DIRTY (Nz=100)','DIRTY (Nz=50)',
                           'DIRTY (Nz=10)','DIRTY (Nz=6)','DIRTY (Nz=3)']
        modnames = ['dirty_stau_0.00250','dirty_stau_0.00500','dirty_stau_0.01000','dirty_stau_0.05000',
                    'dirty_stau_0.10000','dirty_stau_0.25000','dirty_stau_1.00000']
        imodnames = ['dirty_stau/' + modname + '_slab_eff' for modname in modnames]
        scomp = -2
    elif args.nbinz:
        nbinzs = ['10','20','50','100','200','500']
        moddisplaynames = ['DIRTY (Nz='+nbinz+')' for nbinz in reversed(nbinzs)]
        modnames = ['dirty_nbinz_'+nbinz for nbinz in reversed(nbinzs)]
        imodnames = ['dirty_nbinz/' + modname + '_slab_eff' for modname in modnames]
        scomp = -2
    elif args.nphot:
        nphots = ['3.2e5','1e6','3.2e6','1e7','3.2e7','1e8']
        moddisplaynames = ['DIRTY (N='+nphot+')' for nphot in reversed(nphots)]
        modnames = ['dirty_nphot_'+nphot for nphot in reversed(nphots)]
        imodnames = ['dirty_nphot/' + modname + '_slab_eff' for modname in modnames]
        scomp = -2
    elif args.mscat:
        mscats = ['1','5','10','20','50','75','100','150','200','300']
        moddisplaynames = ['DIRTY (mscat=' + mscat + ')' for mscat in reversed(mscats)]
        modnames = ['dirty_mscat_' + mscat for mscat in reversed(mscats)]
        imodnames = ['dirty_mscat/' + modname + '_slab_eff' for modname in modnames]
        scomp = -2
    elif args.econs:
        econtargs = ['1.0','0.32','0.1','0.032','0.01','0.0032','0.001']
        moddisplaynames = ['DIRTY (econs='+econtarg+')' for econtarg in reversed(econtargs)]
        modnames = ['dirty_econs_'+econtarg for econtarg in reversed(econtargs)]
        imodnames = ['dirty_econs/' + modname + '_slab_eff' for modname in modnames]
        scomp = -2
    elif args.nz:
        moddisplaynames = ['SKIRT (Nz=400)','SKIRT (Nz=200)','SKIRT (Nz=100)','SKIRT (Nz=30)','SKIRT (Nz=10)','SKIRT (Nz=5)']
        modnames = ['skirtnz400','skirtnz200','skirtnz100','skirtnz030','skirtnz010','skirtnz005']
        imodnames = ['skirtnz/' + modname + '_slab_eff' for modname in modnames]
        scomp = -2
    elif args.wr:
        moddisplaynames = ['SKIRT (wr=1e8)','SKIRT (wr=1e7)','SKIRT (wr=1e6)','SKIRT (wr=1e5)','SKIRT (wr=1e4)','SKIRT (wr=1e3)']
        modnames = ['skirtwr1e8','skirtwr1e7','skirtwr1e6','skirtwr1e5','skirtwr1e4','skirtwr1e3']
        imodnames = ['skirtwr/' + modname + '_slab_eff' for modname in modnames]
        scomp = -2
    else:
        moddisplaynames = ['CRT','DART-ray','DIRTY','Hyperion','SKIRT','SOC','TRADING']
        modnames = ['crt','dartr','dirty','hyper','skirt','SOC','tradi']
        imodnames = [modname + '/' + modname + '_slab_eff' for modname in modnames]
        scomp = -2

    for angle in angles:
        for tau in taus:
            if args.image:
                for wave in waves:
                    plot_imagegrid(imodnames, moddisplaynames, wave, tau, angle,
                                   save_eps=args.eps, save_png=args.png, comp_index=scomp)
            else:
                plot_decompose_sed(imodnames, moddisplaynames, tau, angle, save_eps=args.eps, save_png=args.png,
                                   single_comp=scomp, max_plot_diff=mplot_diff)
