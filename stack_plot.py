#!/usr/bin/env python

import re
import csv
import os
import sys
import uproot
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import mplhep as hep
import argparse

from load_data import load_data
from pprint import pprint


pjoin = os.path.join

# Pretty labels for legend for each process
pretty_labels = {
    'ZJetsToNuNu' : r'QCD $Z\rightarrow \nu \nu$',
    'EWKZNuNu'    : r'EWK $Z\rightarrow \nu \nu$',
    'DYJetsToLL'  : r'QCD $Z\rightarrow \ell \ell$',
    'EWKZLL'      : r'EWK $Z\rightarrow \ell \ell$',
    'WJetsToLNu'  : r'QCD $W\rightarrow \ell \nu$',
    'EWKW'        : r'EWK $W\rightarrow \ell \nu$'
}

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--region', help='The region to be plotted.')
    parser.add_argument('--noCuts', help='Plot without any additional cuts applied.', action='store_true')
    args = parser.parse_args()
    return args

def stack_plot(inpath, outtag, process_list, csv_file, selection_dicts, region, variable='mjj', include_qcd_estimation=False, qcd_estimation=None):
    '''
    Create a stack plot for the processes specified.
    ==================
    ARGUMENTS:
    ==================
    inpath                 : The path containing input ROOT files
    outtag                 : Output tag to name the output directory 
    process_list           : List of physics processes to be plotted
    csv_file               : The CSV file containing XS + sumw information for each dataset
    selection_dicts        : List of dictionaries, each containing information about a selection.
    region                 : Region to be plotted (A,B,C or D for ABCD method).
    variable               : The variable of interest, by defualt it is mjj.
    include_qcd_estimation : If set to True, include the QCD estimation in the MC stack (False by default).
                             One must provide the QCD estimation as an array if this flag is set to True. 
    qcd_estimation         : If include_qcd_estimation is set to True, provide the QCD estimation as an array (None by default).
    '''
    # Check about the QCD estimation
    if include_qcd_estimation and (qcd_estimation is None):
        raise RuntimeError('Please specify the QCD estimation for plotting.')
    
    print(f'MSG% Starting job, region: {region}')
    # Obtain the histograms for each process specified
    histograms = {}
    data = {}
    # Add the QCD estimation to the MC stack if requested
    if include_qcd_estimation:
        histograms['QCD'] = qcd_estimation

    for process in process_list:
        print(f'MSG% Obtaining histogram for {process}')
        h, bins = load_data(inpath, process, csv_file, variable, selection_dicts)
        if process != 'MET':
            histograms[process] = h
        else:
            data[process] = h

    # Construct the stack plot from the processes
    stacked_histos = np.vstack(histograms.values())
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    # Plot MC
    labels = []
    if include_qcd_estimation:
        labels.append('QCD')

    for process in process_list:
        if process == 'MET':
            continue
        if process in pretty_labels.keys():
            labels.append(pretty_labels[process])
        else:
            labels.append(process)

    hep.histplot(stacked_histos, bins, ax=ax, stack=True, label=labels, binwnorm=True, histtype='fill')
    # Plot data
    kwargs = {
        'linestyle' : 'none',
        'markersize' : 10.,
        'color' : 'k'
    }
    data = np.vstack(data.values())
    hep.histplot(data, bins, ax=ax, label='Data', binwnorm=True, histtype='errorbar', **kwargs)
    ax.legend()
    ax.set_ylabel('Events / Bin Width')
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 1e5)    

    # Plot the data/MC ratio
    total_mc = np.sum(stacked_histos, axis=0)
    ratio = data / total_mc
    hep.histplot(ratio, bins, ax=rax, histtype='errorbar', **kwargs)

    rax.set_ylabel('Data / MC')
    rax.set_ylim(0.5,1.5)
    rax.set_xlabel(r'$M_{jj} \ (GeV)$')

    loc1 = matplotlib.ticker.MultipleLocator(base=0.2)
    loc2 = matplotlib.ticker.MultipleLocator(base=0.1)
    rax.yaxis.set_major_locator(loc1)
    rax.yaxis.set_minor_locator(loc2)
    rax.grid(axis='y',which='minor',linestyle='--')
    rax.grid(axis='y',which='major',linestyle='--')

    # Save figure
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    qcd_estimation_suffix = '_withQCD' if include_qcd_estimation else ''

    if region in ['A', 'B', 'C', 'D']:
        outfile = f'stack_plot_region{region}{qcd_estimation_suffix}.pdf'
    else:
        outfile = f'stack_plot_{region}{qcd_estimation_suffix}.pdf'

    outpath = pjoin(outdir, outfile)
    print(f'MSG% File saved: {outpath}')
    fig.savefig(outpath)

    # Calculate the excess data for each mjj bin
    excess_data = (data - total_mc)[0]
    return excess_data, bins

def main():
    # Path to ROOT files
    # inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees'
    inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/2020-06-30_nodphijj'
    csv_file = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/csv/xs_sumw.csv'

    outtag = os.path.basename(inpath)

    # List of processes to be plotted
    process_list = ['DYJetsToLL', 'Top', 'Diboson', 'EWKW', 'EWKZLL', 'EWKZNuNu', 'WJetsToLNu', 'ZJetsToNuNu', 'MET']

    # Region definitions for ABCD method
    # Region A: dphijj > 1.5 & 1.0 < dPhiTrailingJetMet < 2.3
    # Region B: dphijj > 1.5 & dPhiTrailingJetMet > 2.3
    # Region C: dphijj < 1.5 & 1.0 < dPhiTrailingJetMet < 2.3
    # Region D: dphijj < 1.5 & dPhiTrailingJetMet > 2.3

    args = parse_cli()

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
        ],
        # Cuts for signal region
        'signal': [
            {'variable' : 'dphijj', 'low' : None, 'high' : 1.5}
        ]
    }

    # Pick the relevant selection for the region being specified (if specified)
    if args.region:
        region = args.region
        if region in ['A', 'B', 'C', 'D']:
            selection_dicts = selections_by_region[f'region {region}']
        else:
            selection_dicts = selections_by_region[region] 
    elif (not args.region and args.noCuts):
        region = 'noCuts'
        selection_dicts = None
    else:
        raise RuntimeError('Either specify a region via --region option or specify --noCuts.')

    excess_data, bins = stack_plot(inpath, outtag, process_list, csv_file, selection_dicts=selection_dicts, region=region)

if __name__ == '__main__':
    main()

