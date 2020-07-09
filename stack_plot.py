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
from helper_classes import Style, Selection
from pprint import pprint

pjoin = os.path.join

# Load in the classes holding information about the plots
sty = Style()
xlabels       = sty.xlabels
fig_titles    = sty.fig_titles
pretty_labels = sty.pretty_labels

# Set the selection variables and thresholds
selection_vars = ['dphijj', 'max(neEmEF)']
thresholds = [1.5, 0.8]
sel = Selection(variables=selection_vars, thresholds=thresholds)

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', help='The tree version to be used as inputs, defualt is 09Jul20.', default='09Jul20')
    parser.add_argument('--variables', help='The list of variables to be plotted, default is mjj.', nargs='*', default='mjj')
    parser.add_argument('--region', help='The region to be plotted.')
    parser.add_argument('--noCuts', help='Plot without any additional cuts applied.', action='store_true')
    parser.add_argument('--include_qcd_mc', help='Include the QCD MC in the stack plot.', action='store_true')
    args = parser.parse_args()
    return args

def stack_plot(inpath, outtag, process_list, csv_file, selection_dicts, 
        region, variable='mjj', include_qcd_estimation=False,  
        qcd_estimation=None, include_qcd_mc=False
        ):
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
    include_qcd_mc         : If set to True, QCD MC will be included in the stack plot.
    '''
    # Check about the QCD estimation
    if include_qcd_estimation and (qcd_estimation is None):
        raise RuntimeError('Please specify the QCD estimation for plotting.')
    
    print(f'MSG% Starting job -- Region: {region}, Variable: {variable}')
    # Obtain the histograms for each process specified
    histograms = {}
    data = {}
    # Add the QCD estimation to the MC stack if requested
    if include_qcd_estimation:
        histograms['QCD'] = qcd_estimation

    # If we want to include the QCD MC in the stack, add it to the beginning of the process list
    if include_qcd_mc:
        process_list.insert(0,'QCD')

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

    if region in ['A', 'B', 'C', 'D']:
        ax.set_title(fig_titles[f'region {region}'])
    else:
        ax.set_title(fig_titles[region])

    # Aesthetics: Put edge colors
    handles, labels = ax.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        if label == 'Data':
            continue
        handle.set_linestyle('-')
        handle.set_edgecolor('k')

    # Plot the data/MC ratio
    total_mc = np.sum(stacked_histos, axis=0)
    ratio = data / total_mc
    hep.histplot(ratio, bins, ax=rax, histtype='errorbar', **kwargs)

    rax.set_ylabel('Data / MC')
    rax.set_ylim(0.5,1.5)
    rax.set_xlabel(xlabels[variable])

    loc1 = matplotlib.ticker.MultipleLocator(base=0.2)
    loc2 = matplotlib.ticker.MultipleLocator(base=0.1)
    rax.yaxis.set_major_locator(loc1)
    rax.yaxis.set_minor_locator(loc2)
    rax.grid(axis='y',which='minor',linestyle='--')
    rax.grid(axis='y',which='major',linestyle='--')

    xlim = rax.get_xlim()
    rax.plot(xlim, [1., 1.], 'r--')
    rax.set_xlim(xlim)

    # Save figure
    if region in ['A', 'B', 'C', 'D'] or include_qcd_estimation:
        selection_tag = sel.selection_tag
        outdir = f'./output/{outtag}/qcd_estimation/{selection_tag}'
    elif include_qcd_mc:
        outdir = f'./output/{outtag}/with_qcd_mc'
    else:
        outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    qcd_estimation_suffix = '_withQCD' if include_qcd_estimation else ''

    if region in ['A', 'B', 'C', 'D']:
        outfile = f'stack_plot_region{region}_{variable}{qcd_estimation_suffix}.pdf'
    else:
        outfile = f'stack_plot_{region}_{variable}{qcd_estimation_suffix}.pdf'

    outpath = pjoin(outdir, outfile)
    print(f'MSG% File saved: {outpath}')
    fig.savefig(outpath)

    # Calculate the excess data for each mjj bin
    excess_data = (data - total_mc)[0]
    return excess_data, bins

def main():
    args = parse_cli()
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

    # List of processes to be plotted
    process_list = ['DYJetsToLL', 'Top', 'Diboson', 'EWKW', 'EWKZLL', 'EWKZNuNu', 'WJetsToLNu', 'ZJetsToNuNu', 'MET']

    # Pick the relevant selection for the region being specified (if specified)
    if args.region:
        region = args.region
        if region in ['A', 'B', 'C', 'D']:
            selection_dicts = sel.selections_by_region[f'region {region}']
        else:
            selection_dicts = sel.selections_by_region[region]
    elif (not args.region and args.noCuts):
        region = 'noCuts'
        selection_dicts = None
    else:
        raise RuntimeError('Either specify a region via --region option or specify --noCuts.')

    # If a variable list is given, loop over each variable and make plots
    if isinstance(args.variables, list): 
        for variable in args.variables:
            excess_data, bins = stack_plot(inpath, outtag, process_list, csv_file, 
                                variable=variable,
                                selection_dicts=selection_dicts, 
                                region=region,
                                include_qcd_mc=args.include_qcd_mc
                                )
    # The case where only a single variable is specified
    else:
        variable = args.variables
        excess_data, bins = stack_plot(inpath, outtag, process_list, csv_file, 
                            variable=variable,
                            selection_dicts=selection_dicts, 
                            region=region,
                            include_qcd_mc=args.include_qcd_mc
                            )

if __name__ == '__main__':
    main()

