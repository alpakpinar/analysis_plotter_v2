#!/usr/bin/env python

import re
import csv
import os
import sys
from pprint import pprint
import uproot
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import mplhep as hep
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

def stack_plot(inpath, outtag, process_list, csv_file, selection_dicts, variable='mjj'):
    '''
    Create a stack plot for the processes specified.
    ==================
    ARGUMENTS:
    ==================
    inpath          : The path containing input ROOT files
    outtag          : Output tag to name the output directory 
    process_list    : List of physics processes to be plotted
    csv_file        : The CSV file containing XS + sumw information for each dataset
    selection_dicts : List of dictionaries, each containing information about a selection.
    variable        : The variable of interest, by defualt it is mjj.
    '''
    # Obtain the histograms for each process specified
    histograms = {}
    data = {}
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
    
    cut_label_for_file = ''
    if selection_dicts is not None:
        for selection_dict in selection_dicts:
            cut_label_for_file += selection_dict['file_label'] 
    else:
        cut_label_for_file = '_sr_cuts'

    outfile = f'stack_plot{cut_label_for_file}.pdf'
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
    print(f'MSG% CSV file containing XS and sumw: {csv_file}')

    outtag = os.path.basename(inpath)

    # List of processes to be plotted
    process_list = ['DYJetsToLL', 'Top', 'Diboson', 'EWKW', 'EWKZLL', 'EWKZNuNu', 'WJetsToLNu', 'ZJetsToNuNu', 'MET']

    selection_dicts = [
        # {'variable' : 'dphijj', 'low' : None, 'high' : 1.5, 'label' : 'dphijj<1.5', 'file_label' : '_dphijj_smallerThan_1_5'},
        # {'variable' : 'dphijj', 'low' : 1.5, 'high' : None, 'label' : 'dphijj>1.5', 'file_label' : '_dphijj_largerThan_1_5'},
        # {'variable' : 'dPhiTrailingJetMet', 'low' : 1.0, 'high' : 2.3, 'label' : '1.0<dPhiTrailingJetMet<2.3', 'file_label' : '_dPhiTrailingJetMet_between_1_0_and_2_3'}
        {'variable' : 'dPhiTrailingJetMet', 'low' : 2.3, 'high' : None, 'label' : 'dPhiTrailingJetMet>2.3', 'file_label' : '_dPhiTrailingJetMet_largerThan_2_3'}
    ]

    excess_data, bins = stack_plot(inpath, outtag, process_list, csv_file, selection_dicts=selection_dicts)

if __name__ == '__main__':
    main()

