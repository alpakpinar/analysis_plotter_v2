#!/usr/bin/env python

# =============================
# Script to get the ratios of excess data events in different 
# regions as a function of mjj, for QCD estimation studies. 
# =============================

import os
import sys
import re
import numpy as np
import mplhep as hep
import argparse
import uproot
from matplotlib import pyplot as plt
from pprint import pprint
from stack_plot import stack_plot

pjoin = os.path.join

selections_by_region = {
    'region A' : [
        {'variable' : 'dphijj', 'low' : 1.5, 'high' : None},
        {'variable' : 'dPhiTrailingJetMet', 'low' : 1.0, 'high' : 2.3}
    ],
    'region B' : [
        {'variable' : 'dphijj', 'low' : 1.5, 'high' : None},
        {'variable' : 'dPhiTrailingJetMet', 'low' : 2.3, 'high' : None}
    ],
    'region C' : [
        {'variable' : 'dphijj', 'low' : None, 'high' : 1.5},
        {'variable' : 'dPhiTrailingJetMet', 'low' : 1.0, 'high' : 2.3}
    ],
    'region D' : [
        {'variable' : 'dphijj', 'low' : None, 'high' : 1.5},
        {'variable' : 'dPhiTrailingJetMet', 'low' : 2.3, 'high' : None}
    ]
}

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--region1', help='First region to be run over.')
    parser.add_argument('--region2', help='Second region to be run over.')
    args = parser.parse_args()
    return args

def get_ratio_of_excess_data(inpath, outtag, region1, region2, process_list, csv_file, variable='mjj', save_to_root=True):    
    '''Get the ratio of excess data events (over MC) in two regions, region1 and region2.'''
    # Call the stack_plot function to get the excess events in each region
    excess_events = {}
    for region in [region1, region2]:
        selection_dicts = selections_by_region[f'region {region}'] if region != 'noCuts' else None
        excess_events[region], bins = stack_plot(inpath, outtag, process_list, 
                                                 csv_file, 
                                                 selection_dicts=selection_dicts,
                                                 region=region
                                                 )

        # If excess events < 0, can set them to zero since we are not interested with those bins
        excess_events[region][excess_events[region] < 0] = 0.

    # Plot the excess events for each region as a function of mjj
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    hep.histplot(excess_events[region1], bins, ax=ax, label=f'Region {region1}', histtype='step')
    hep.histplot(excess_events[region2], bins, ax=ax, label=f'Region {region2}', histtype='step')

    ax.set_ylabel('Excess number of events')
    ax.legend(title='Regions')

    # Calculate and plot the ratio: region1/region2
    kwargs = {
        'linestyle' : 'none',
        'markersize' : 10.,
        'color' : 'k'
    }
    ratio = excess_events[region1] / excess_events[region2]
    hep.histplot(ratio, bins, ax=rax, histtype='errorbar', **kwargs)

    rax.grid(True)
    rax.set_xlabel(r'$M_{jj} \ (GeV)$')
    rax.set_ylabel(f'{region1} / {region2}')

    # Save the figure
    outdir = f'./output/{outtag}/excess_events'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = f'excess_events_regions_{region1}_{region2}.pdf'

    outpath = pjoin(outdir, outfile)
    fig.savefig(outpath)
    print(f'MSG% Figure saved: {outpath}')

    # Save ratio to a ROOT file if requested (default behavior)
    if save_to_root:
        # If the ROOT file is already there do not create a new one
        root_filepath = pjoin(outdir, 'excess_event_ratios.root')
        if os.path.isfile(root_filepath):
            f = uproot.open(root_filepath)
            print(f'MSG% Modifying the already existing ROOT file: {root_filepath}')
        # If the ROOT file doesn't exist, create a new one
        else:
            f = uproot.recreate(root_filepath)
            print(f'MSG% ROOT file created at: {root_filepath}')

        histo_tag = f'ratio_regions_{region1}_{region2}'
        f[histo_tag] = (ratio, bins)        
        print(f'MSG% Histogram saved to ROOT file: {histo_tag}')

    # Return the ratio and the corresponding binning
    return ratio, bins

def main():
    args = parse_cli()

    # Warn the user if regions are not specified
    if (not args.region1) or (not args.region2):
        raise RuntimeError('Please specify the two regions for the ratio via --region1 and --region2 options.')

    # Path to ROOT files
    # inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees'
    inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/2020-06-30_nodphijj'
    csv_file = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/csv/xs_sumw.csv'

    outtag = os.path.basename(inpath)

    # List of processes to be plotted
    process_list = ['DYJetsToLL', 'Top', 'Diboson', 'EWKW', 'EWKZLL', 'EWKZNuNu', 'WJetsToLNu', 'ZJetsToNuNu', 'MET']

    excess_ratio, bins = get_ratio_of_excess_data(inpath, outtag, region1=args.region1, region2=args.region2, process_list=process_list, csv_file=csv_file)

if __name__ == '__main__':
    main()
