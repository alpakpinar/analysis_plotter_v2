#!/usr/bin/env python

# ==========================
# Script to plot VBF signal as a function of mjj (the fit variable)
# before and after some cuts are applied, and understand the effect
# of the cut on the signal yield in the mjj spectrum.
# ==========================

import re
import csv
import os
import sys
import uproot
import numpy as np
from matplotlib import pyplot as plt
import mplhep as hep
import argparse

from load_data import load_data
from helper_classes import Cut
from pprint import pprint

pjoin = os.path.join

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', help='The tree version to be used as inputs, defualt is 09Jul20.', default='09Jul20')
    parser.add_argument('--cuts', help='List of cuts to be applied to the signal.', nargs='*', default=[])
    args = parser.parse_args()
    return args

def calc_percent_loss(h_nom, h_all, bins):
    '''Given the old and new histograms, calculate the percentage loss in yield.'''
    bin_widths = np.diff(bins)
    total_yield_nom = np.sum(h_nom*bin_widths)
    total_yield_all = np.sum(h_all*bin_widths)

    total_percent_loss = ( np.abs(total_yield_nom - total_yield_all) / total_yield_nom ) * 100

    # Calculate the loss of yield for mjj > 2 TeV
    yield_large_mjj_nom = np.sum(h_nom[-2:]*bin_widths[-2:])
    yield_large_mjj_all = np.sum(h_all[-2:]*bin_widths[-2:])

    percent_loss_large_mjj = ( np.abs(yield_large_mjj_nom - yield_large_mjj_all) / yield_large_mjj_nom ) * 100

    yield_losses = {
        'total' : total_percent_loss,
        'large_mjj' : percent_loss_large_mjj
    }

    return yield_losses

def plot_signal_comp(inpath, outtag, csv_file, extra_cut_names, variable='mjj'):
    '''Plot a before/after comparison plot for VBF signal MC, given the list of cuts to be applied.'''
    # Nominal selections for signal region
    nom_selection_list = [
        Cut('dphijj', low_thresh=None, high_thresh=1.5),
        Cut('recoil_pt', low_thresh=250, high_thresh=None)
    ]
    
    # Load the data with nominal selections
    h_nom, bins = load_data(inpath, process='VBF', csv_file=csv_file, variable='mjj', cuts=nom_selection_list)

    # Add in the cuts for region D
    nom_selection_list.extend([
        Cut('max(neEmEF)', low_thresh=None, high_thresh=0.7)
    ])

    # List with the additional selections
    extra_selection_list = []

    additional_cuts = {
        'jet_eta'  : Cut('leadak4_trailak4_eta', low_thresh=None, high_thresh=2.5),
        'met_dphi' : Cut('dPhi_TkMET_PFMET', low_thresh=None, high_thresh=1.0),
        'leading_jet_pt' : Cut('leadak4_pt', low_thresh=100, high_thresh=None),
        'leading_jet_pt120' : Cut('leadak4_pt', low_thresh=120, high_thresh=None),
        # Cuts to apply only if none of the two leading jets is in HF
        'met_dphi_noHF'  : Cut('dPhi_TkMET_PFMET', low_thresh=None, high_thresh=0.75, special_apply='noJetInHF'),
        'leading_jet_pt_noHF'  : Cut('leadak4_pt', low_thresh=100, high_thresh=None, special_apply='noJetInHF'),
        'leading_jet_pt120_noHF'  : Cut('leadak4_pt', low_thresh=120, high_thresh=None, special_apply='noJetInHF'),
        # Cuts to apply only if one of the two leading jets is in endcap 
        'met_dphi_jetInEndcap' : Cut('dPhi_TkMET_PFMET', low_thresh=None, high_thresh=0.75, special_apply='oneJetInEndcap'),
        'leading_jet_pt_jetInEndcap' : Cut('leadak4_pt', low_thresh=100, high_thresh=None, special_apply='oneJetInEndcap'),
        'leading_jet_pt120_jetInEndcap' : Cut('leadak4_pt', low_thresh=120, high_thresh=None, special_apply='oneJetInEndcap'),
        # Cuts with endcap coverage slightly increased to |eta| 3.2
        'met_dphi_jetInEndcap_v2' : Cut('dPhi_TkMET_PFMET', low_thresh=None, high_thresh=0.75, special_apply='oneJetInEndcap', change_endcap_def=True),
        'leading_jet_pt_jetInEndcap_v2' : Cut('leadak4_pt', low_thresh=100, high_thresh=None, special_apply='oneJetInEndcap', change_endcap_def=True),
        'leading_jet_pt120_jetInEndcap_v2' : Cut('leadak4_pt', low_thresh=120, high_thresh=None, special_apply='oneJetInEndcap', change_endcap_def=True)
    }
    
    if len(extra_cut_names) != 0:
        for cutname, cut in additional_cuts.items():
            if cutname in extra_cut_names:
                extra_selection_list.append(cut)

    all_selection_list = nom_selection_list + extra_selection_list

    # Load the data with the extra selections applied
    h_all, bins = load_data(inpath, process='VBF', csv_file=csv_file, variable='mjj', cuts=all_selection_list)

    # Plot the comparison between the two 
    fig, ax = plt.subplots()
    hep.histplot(h_nom, bins, ax=ax, label='VBF (Nominal Cuts)', histtype='step')
    hep.histplot(h_all, bins, ax=ax, label='VBF (All Cuts)', histtype='step')

    ax.set_xlabel(r'$M_{jj} \ (GeV)$')
    ax.set_xlim(200,3500)
    ax.set_ylabel('Events / Bin Width')
    ax.set_yscale('log')
    ax.set_ylim(1e1,1e3)
    ax.legend()

    # Calculate the percentage loss in signal yield
    yield_losses = calc_percent_loss(h_nom, h_all, bins)
    total_percent_loss = yield_losses['total']
    percent_loss_large_mjj = yield_losses['large_mjj']

    # Display the losses in the figure
    text1 = f'Total loss: {total_percent_loss:.2f}%'
    text2 = f'Loss for mjj > 2TeV: {percent_loss_large_mjj:.2f}%'

    ax.text(0.6, 0.8, text1, transform=ax.transAxes)
    ax.text(0.6, 0.75, text2, transform=ax.transAxes)

    # Save figure
    outdir = f'./output/{outtag}/vbf_comparison'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if len(extra_cut_names) != 0:
        filename = f'{"_".join(extra_cut_names)}.pdf'    
    else:
        filename = 'noExtraCuts.pdf'
    outpath = pjoin(outdir, filename)
    fig.savefig(outpath)

    print(f'MSG% File saved: {outpath}')

def main():
    args = parse_cli()
    # Get the list of cuts to be applied
    extra_cut_names = args.cuts

    # Path to ROOT files, by default use the latest ones (09Jul20), if specified use 30Jun20 instead.
    if args.version == '09Jul20':
        inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/09Jul20'
    elif args.version == '05Jul20':
        inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/2020-07-05_nodphijj'
    elif args.version == '30Jun20':
        inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/2020-06-30_nodphijj'

    print(f'MSG% Using trees from version: {args.version}')
    # Path to CSV file containing XS + sumw information for every dataset
    csv_file = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/csv/xs_sumw.csv'

    outtag = os.path.basename(inpath)

    plot_signal_comp(inpath, outtag, csv_file, extra_cut_names)

if __name__ == '__main__':
    main()
