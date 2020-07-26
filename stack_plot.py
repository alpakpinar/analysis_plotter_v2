#!/usr/bin/env python

import re
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import mplhep as hep
import argparse

from load_data import load_data
from helper_classes import Style, Selection
from pprint import pprint

pjoin = os.path.join

def load_style_and_selection(additional_cuts, categorization=None):
    '''
    Load classes that contain information about plotting style and selections.
    This function will only be called if this script is the main script being called.
    '''
    # Load in the classes holding information about the plots
    sty = Style()

    # Set the selection variables and thresholds
    selection_vars = ['dphijj', 'max(neEmEF)']
    # selection_vars = ['dphijj', 'dPhi_TkMET_PFMET']
    thresholds = [1.5, 0.7]
    sel = Selection(variables=selection_vars, thresholds=thresholds, apply_cuts=additional_cuts, categorization=categorization)

    return sty, sel

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', help='The tree version to be used as inputs, default is 09Jul20.', default='09Jul20')
    parser.add_argument('--variables', help='The list of variables to be plotted, default is mjj.', nargs='*', default='mjj')
    parser.add_argument('--region', help='The region to be plotted.')
    parser.add_argument('--noCuts', help='Plot without any additional cuts applied.', action='store_true')
    parser.add_argument('--include_qcd_mc', help='Include the QCD MC in the stack plot.', action='store_true')
    parser.add_argument('--additionalCuts', help='Additional cuts to apply on all ABCD regions.', nargs='*', default=['recoil'])
    parser.add_argument('--jesVariation', help='JES variation to look at, can be central, up or down.')
    parser.add_argument('--categorization', help='The categorization to plot (e.g. Trk-EE), by default no categorization cut will be applied.')
    args = parser.parse_args()
    return args

def stack_plot(inpath, outtag, process_list, csv_file, cuts, sty, sel, region, 
        variable='mjj', include_qcd_estimation=False, plot_signal=True, qcd_estimation=None, 
        include_qcd_mc=False, eta_binning='very_fine', output_dir_tag=None, jes_variation='central',
        categorization=None
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
    cuts                   : List of Cut objects, each containing information about a cut.
    sty                    : The Style object to be passed in for plotting.
    sel                    : The Selection object to be passed in for plotting.
    region                 : Region to be plotted (A,B,C or D for ABCD method).
    variable               : The variable of interest, by defualt it is mjj.
    include_qcd_estimation : If set to True, include the QCD estimation in the MC stack (False by default).
                             One must provide the QCD estimation as an array if this flag is set to True. 
    plot_signal            : If set to True, include the signal MC template in the stack plot (True by default)
    qcd_estimation         : If include_qcd_estimation is set to True, provide the QCD estimation as an array (None by default).
    include_qcd_mc         : If set to True, QCD MC will be included in the stack plot.
    eta_binning            : The eta binning to be used for the TF calculation in ABCD method, defaults to "very_fine"
    output_dir_tag         : The tag to be used for the naming of the output directory, depending on the eta binning being used in TF calculation.
    jes_variation          : JES variation to look at. By default it is central (no variation).
    categorization         : The categorization according to the two leading jets (by default, None).
    '''
    # Check about the QCD estimation
    if include_qcd_estimation and (qcd_estimation is None):
        raise RuntimeError('Please specify the QCD estimation for plotting.')
    
    print(f'MSG% Starting job -- Region: {region}, Variable: {variable}')
    if categorization is not None:
        print(f'MSG% Jet categorization: {categorization}')
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
        # Load the data from ROOT files: Get the histograms for each process + binning
        # To smooth out the TF as a function of jet eta in QCD estimation, use coarser eta binning if requested 
        h, bins = load_data(inpath, process, csv_file, variable, cuts, eta_binning=eta_binning, jes_variation=jes_variation)

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
        if process in sty.pretty_labels.keys():
            labels.append(sty.pretty_labels[process])
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
    ax.set_ylabel('Events / Bin Width')
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 1e8)    

    # Include signal in the plot if requested
    if plot_signal:
        signal, bins = load_data(inpath, process='VBF', csv_file=csv_file, variable=variable, cuts=cuts, eta_binning=eta_binning)
        hep.histplot(signal, bins, ax=ax, label='VBF Signal', binwnorm=True, histtype='step')

    if region in ['A', 'B', 'C', 'D']:
        ax.set_title(sel.fig_titles[f'region {region}'])
    else:
        ax.set_title(sel.fig_titles[region])

    # Aesthetics: Put edge colors
    handles, labels = ax.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        if label == 'Data':
            continue
        elif label == 'VBF Signal':
            handle.set_color('k')
            handle.set_linewidth(2)
            continue
        handle.set_linestyle('-')
        handle.set_edgecolor('k')

    ax.legend(ncol=2)

    # Plot the data/MC ratio
    total_mc = np.sum(stacked_histos, axis=0)
    ratio = data / total_mc
    hep.histplot(ratio, bins, ax=rax, histtype='errorbar', **kwargs)

    rax.set_ylabel('Data / MC')
    rax.set_ylim(0.5,1.5)
    rax.set_xlabel(sty.xlabels[variable])

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
        additional_selection_tag = sel.additional_selection_tag  
        if output_dir_tag:
            outdir = f'./output/{outtag}/qcd_estimation/{selection_tag}/{additional_selection_tag}/{output_dir_tag}'
        else:
            outdir = f'./output/{outtag}/qcd_estimation/{selection_tag}/{additional_selection_tag}'
    elif region =='signal' and jes_variation is not None:
        outdir = f'./output/{outtag}/{jes_variation}'
    elif include_qcd_mc:
        outdir = f'./output/{outtag}/with_qcd_mc'
    else:
        outdir = f'./output/{outtag}'
    # Add the categorized plot in a separate sub-directory
    if categorization is not None:
        additional_selection_tag = sel.additional_selection_tag
        if additional_selection_tag != 'recoil':
            outdir += f'/categorized/{categorization}/{additional_selection_tag}'  
        else:
            outdir += f'/categorized/{categorization}'  
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
    additional_cuts = args.additionalCuts
    sty, sel = load_style_and_selection(additional_cuts, categorization=args.categorization)

    # Path to ROOT files, by default use the latest ones (09Jul20), if specified use 30Jun20 instead.
    if args.jesVariation is not None or args.version == '25Jul20':
        # Use the trees with JES variations recorded
        inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/25Jul20_JES'
    else:
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
            cuts = sel.selections_by_region[f'region {region}']
        else:
            cuts = sel.selections_by_region[region]
    elif (not args.region and args.noCuts):
        region = 'noCuts'
        cuts = None
    else:
        raise RuntimeError('Either specify a region via --region option or specify --noCuts.')

    # If a variable list is given, loop over each variable and make plots
    if isinstance(args.variables, list): 
        for variable in args.variables:
            excess_data, bins = stack_plot(inpath, outtag, process_list, csv_file, 
                                variable=variable,
                                sty=sty, sel=sel,
                                cuts=cuts, 
                                region=region,
                                include_qcd_mc=args.include_qcd_mc,
                                jes_variation=args.jesVariation,
                                categorization=args.categorization
                                )
    # The case where only a single variable is specified
    else:
        variable = args.variables
        excess_data, bins = stack_plot(inpath, outtag, process_list, csv_file, 
                            variable=variable,
                            sty=sty, sel=sel,
                            cuts=cuts, 
                            region=region,
                            include_qcd_mc=args.include_qcd_mc,
                            jes_variation=args.jesVariation,
                            categorization=args.categorization
                            )

if __name__ == '__main__':
    main()

