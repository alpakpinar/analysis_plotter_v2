#!/usr/bin/env python

# =============================
# Script to get the QCD estimate values by using ABCD method. 
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
from helper_classes import Style, Selection

pjoin = os.path.join

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', help='The tree version to be used as inputs, defualt is 29Jul20.', default='29Jul20')
    parser.add_argument('--variable', help='The variable for the plotting of QCD template.', default='mjj')
    parser.add_argument('--eta_binning', help='The eta binning for the calculation of TF: C/B, can be fine or coarse. By default, very_fine is used.', default='very_fine')
    parser.add_argument('--selection_vars', help='The two variables used to define the ABCD regions.', nargs=2, default=['dphijj', 'max(neEmEF)'])
    parser.add_argument('--thresholds', help='Thresholds on the two variables that are used to define ABCD regions.', nargs=2, type=float, default=[1.5, 0.7])
    parser.add_argument('--additionalCuts', help='Additional cuts to apply on all ABCD regions.', nargs='*', default=['recoil'])
    parser.add_argument('--applyCleaningCuts', help='Apply cleaning cuts on all ABCD regions.', action='store_true')
    parser.add_argument('--veto_hfhf', help='Veto HF-HF events.', action='store_true')
    args = parser.parse_args()
    return args

def get_ratio_of_excess_data(inpath, outtag, region1, region2, process_list, csv_file, sty, sel, variable='mjj', 
            save_to_root=False, eta_binning='very_fine', output_dir_tag=None, apply_cleaning_cuts=False, veto_hfhf=False
            ):    
    '''Get the ratio of excess data events (over MC) in two regions, region1 and region2.'''
    # Call the stack_plot function to get the excess events in each region
    excess_events = {}
    for region in [region1, region2]:
        cuts = sel.selections_by_region[f'region {region}'] if region != 'noCuts' else None
        excess_events[region], bins = stack_plot(inpath, outtag, process_list, 
                                                 csv_file,
                                                 variable=variable,
                                                 sel=sel, sty=sty, 
                                                 cuts=cuts,
                                                 region=region,
                                                 eta_binning=eta_binning,
                                                 output_dir_tag=output_dir_tag,
                                                 apply_cleaning_cuts=apply_cleaning_cuts,
                                                 veto_hfhf=veto_hfhf
                                                 )

        # If excess events < 0, can set them to zero since we are not interested with those bins
        excess_events[region][excess_events[region] < 0] = 0.


    # Plot the excess events for each region as a function of the requested variable
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    hep.histplot(excess_events[region1], bins, ax=ax, label=f'Region {region1}', histtype='step')
    hep.histplot(excess_events[region2], bins, ax=ax, label=f'Region {region2}', histtype='step')

    ax.set_ylabel('Excess number of events')
    ax.set_yscale('log')
    ax.set_ylim(1e-1, 1e5)
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
    rax.set_xlabel(sty.xlabels[variable])
    rax.set_ylabel(f'{region1} / {region2}')
    rax.set_ylim(0,5)

    # Save the figure
    if output_dir_tag:
        outdir = f'./output/{outtag}/qcd_estimation/{sel.selection_tag}/{sel.additional_selection_tag}/{output_dir_tag}'
    else:
        outdir = f'./output/{outtag}/qcd_estimation/{sel.selection_tag}/{sel.additional_selection_tag}'
    if apply_cleaning_cuts:
        outdir += '/with_cleaning_cuts'
    if veto_hfhf:
        outdir += '/veto_hfhf'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = f'excess_events_regions_{region1}_{region2}_{variable}.pdf'

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

def get_qcd_estimate(inpath, outtag, process_list, csv_file, sel, sty, variable='mjj', 
                save_to_root=False, eta_binning='very_fine', output_dir_tag=None, apply_cleaning_cuts=False,
                veto_hfhf=False
                ):
    '''
    Using the ratios between several regions, get the QCD estimate for the signal region.
    ==================
    ARGUMENTS:
    ==================
    inpath           : The path containing input ROOT files
    outtag           : Output tag to name the output directory 
    process_list     : List of physics processes to be plotted 
    csv_file         : The CSV file containing XS + sumw information for each dataset
    sel              : The object containing information on selection.
    sty              : The object containing information on plotting style & labels.
    variable         : The variable of interest, by defualt it is mjj
    save_to_root     : If set to True, save the results into an output ROOT file
    eta_binning      : The eta binning to be used for the TF calculation in ABCD method, defaults to "very_fine"
    output_dir_tag   : The tag to be used for output directory naming, has information about the eta binning used in C/B calculation.
    apply_cleaning_cuts   : Apply cleaning cuts (VecB, VecDPhi) on top of all ABCD regions, by default, this will not be done.
    veto_hfhf        : Veto HF-HF events in the calculation, by default this is not done.
    '''
    # Here, the QCD estimation is calculated as: (C/B) * A 
    # First, get the ratio of C/B, use coarser eta binning for the ratio (to smooth the TF) if requested
    ratio_C_B, bins = get_ratio_of_excess_data(inpath, outtag, region1='C', region2='B', 
                                    variable=variable, 
                                    process_list=process_list, 
                                    csv_file=csv_file,
                                    sty=sty, sel=sel, 
                                    save_to_root=False,
                                    eta_binning=eta_binning,
                                    output_dir_tag=output_dir_tag,
                                    apply_cleaning_cuts=apply_cleaning_cuts,
                                    veto_hfhf=veto_hfhf
                                    )
    # If coarser eta binning is used, resize the ratio array by repetition so that its compatible 
    # to use in arithmetic operations with the other histograms
    if eta_binning == 'fine':
        ratio_C_B = np.repeat(ratio_C_B,5)
    elif eta_binning == 'coarse':
        if 'eta' in variable:
            ratio_C_B = np.repeat(ratio_C_B,10)
    elif eta_binning == 'coarse_largeEta':
        if not 'absEta' in variable:
            raise RuntimeError('This binning is currently only defined for absEta variable!')
        # Extend the last TF to the largest eta bins
        num_bins_to_extend = 8
        ratio_C_B = list(ratio_C_B[:-1]) + list(np.repeat( ratio_C_B[-1], num_bins_to_extend ))

    # Get the excess data events for region A
    excess_events_A, bins = stack_plot(inpath, outtag, 
                                    variable=variable,
                                    process_list=process_list,
                                    sel=sel, sty=sty, 
                                    csv_file=csv_file, 
                                    cuts=sel.selections_by_region['region A'], 
                                    region='A',
                                    output_dir_tag=output_dir_tag,
                                    apply_cleaning_cuts=apply_cleaning_cuts,
                                    veto_hfhf=veto_hfhf
                                    ) 

    # If excess events are smaller than 0, just set them to 0 since we're not interested in those
    excess_events_A[excess_events_A < 0] = 0.

    # Plot the excess events in region A and save
    fig, ax = plt.subplots()
    hep.histplot(excess_events_A, bins, ax=ax, histtype='step')

    ax.set_ylabel('Excess Number of Events')
    ax.set_xlabel(sty.xlabels[variable])
    ax.set_title('Region A')

    ax.set_yscale('log')
    ax.set_ylim(1e-1, 1e4)

    if output_dir_tag:
        outdir = f'./output/{outtag}/qcd_estimation/{sel.selection_tag}/{sel.additional_selection_tag}/{output_dir_tag}'
    else:
        outdir = f'./output/{outtag}/qcd_estimation/{sel.selection_tag}/{sel.additional_selection_tag}'
    
    if apply_cleaning_cuts:
        outdir += '/with_cleaning_cuts'
    if veto_hfhf:
        outdir += '/veto_hfhf'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outfile = f'excess_events_regionA_{variable}.pdf'
    outpath = pjoin(outdir, outfile)
    fig.savefig(outpath)
    print(f'MSG% File saved: {outpath}')
    plt.close(fig)

    # Get the QCD estimate for region D (region of interest)
    bad_value = np.isnan(ratio_C_B) | np.isinf(ratio_C_B)
    ratio_C_B = np.where(bad_value, 1, ratio_C_B)
    qcd_estimation = ratio_C_B * excess_events_A

    # Plot the QCD estimation as a function of mjj
    fig, ax = plt.subplots()
    hep.histplot(qcd_estimation, bins, ax=ax, histtype='step')
    ax.set_xlabel(sty.xlabels[variable])
    ax.set_ylabel('Events')
    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e5)
    ax.set_title('QCD Estimation')
    
    # Save figure
    outpath = pjoin(outdir, f'qcd_estimation_{variable}.pdf')
    fig.savefig(outpath)
    print(f'MSG% File saved: {outpath}')

    # Return the QCD estimations and the corresponding binning
    return qcd_estimation, bins

def stack_plot_with_qcd_estimation(inpath, outtag, variable, sel, sty, process_list, csv_file, qcd_estimation, output_dir_tag=None, apply_cleaning_cuts=False, veto_hfhf=False):
    '''
    Create a stack plot for the signal region with the QCD estimation included.
    Specify the pre-calculated QCD-estimation in qcd_estimation parameter as an array.
    '''
    # Call the stack_plot function with the QCD estimate included
    region = 'D'
    cuts = sel.selections_by_region[f'region {region}']

    stack_plot(inpath, outtag, process_list, csv_file,
               variable=variable, 
               cuts=cuts, 
               sel=sel, sty=sty,
               region=region, 
               include_qcd_estimation=True, 
               qcd_estimation=qcd_estimation,
               output_dir_tag=output_dir_tag,
               apply_cleaning_cuts=apply_cleaning_cuts,
               veto_hfhf=veto_hfhf
               )

def main():
    args = parse_cli()

    # Path to ROOT files, by default use the latest ones (09Jul20), if specified use 05Jul20 or 30Jun20 instead.
    version = args.version
    if version == '29Jul20':
        inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/29Jul20'
    elif version == '09Jul20':
        inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/09Jul20'
    elif version == '05Jul20':
        inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/2020-07-05_nodphijj'
    elif version == '30Jun20':
        inpath = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/2020-06-30_nodphijj'

    print(f'MSG% Using trees from version: {version}')
    variable = args.variable
    print(f'MSG% Getting templates as a function of: {variable}')
    
    # Path to CSV file containing XS + sumw information for every dataset
    csv_file = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/csv/xs_sumw.csv'

    outtag = os.path.basename(inpath)

    # Read variables and thresholds for ABCD method from the command line
    selection_vars = args.selection_vars
    thresholds = args.thresholds

    print(f'MSG% Variables used for ABCD method: {", ".join(selection_vars)}')
    print(f'MSG% Thresholds used for ABCD method: {", ".join(map(str, thresholds))}')

    # Load in the classes holding information about the plots:
    sty = Style()
    
    sel = Selection(variables=selection_vars, thresholds=thresholds, apply_cuts=args.additionalCuts)

    # List of processes to be plotted
    process_list = ['DYJetsToLL', 'Top', 'Diboson', 'EWKW', 'EWKZLL', 'EWKZNuNu', 'WJetsToLNu', 'ZJetsToNuNu', 'MET']

    # Eta binning to be used in transfer factor calculation
    eta_binning = args.eta_binning
    # Plots will go into this output dir
    output_dir_tag = f'tf_eta_binning_{eta_binning}'

    qcd_estimation, bins = get_qcd_estimate(inpath, outtag, 
                                process_list=process_list, 
                                csv_file=csv_file,
                                sel=sel, sty=sty,
                                variable=variable,
                                eta_binning=eta_binning,
                                output_dir_tag=output_dir_tag,
                                apply_cleaning_cuts=args.applyCleaningCuts,
                                veto_hfhf=args.veto_hfhf
                                )

    # Create a stack plot with QCD estimation included
    stack_plot_with_qcd_estimation(inpath, outtag, 
                                process_list=process_list, 
                                csv_file=csv_file,
                                sel=sel, sty=sty, 
                                variable=variable,
                                qcd_estimation=qcd_estimation,
                                output_dir_tag=output_dir_tag,
                                apply_cleaning_cuts=args.applyCleaningCuts,
                                veto_hfhf=args.veto_hfhf
                                )

if __name__ == '__main__':
    main()
