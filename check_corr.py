#!/usr/bin/env python

import os
import sys
from pprint import pprint
import uproot
import numpy as np
from matplotlib import pyplot as plt

pjoin = os.path.join

binnings = {
    'dphijj' : np.linspace(0,3,16),
    'max(neEmEF)' : np.linspace(0,1,11),
}

def plot_2d_histogram(input_root_file, proc, tree_version):
    '''Given the input ROOT file, plot a 2D histogram of dphijj and jet fractions.'''
    events = uproot.open(input_root_file)['sr_vbf']

    # Define signal region mask
    mask = events['met_pt'].array() > 250 
    dphijj_array = events['dphijj'].array()[mask]

    leadak4_neEmEF = events['leadak4_neEmEF'].array()[mask]
    trailak4_neEmEF = events['trailak4_neEmEF'].array()[mask]
    max_neEmEF_array = np.maximum(leadak4_neEmEF, trailak4_neEmEF)

    # Create a 2D histogram
    fig, ax = plt.subplots()
    h, xedges, yedges, im = ax.hist2d(dphijj_array, max_neEmEF_array, bins=[ binnings['dphijj'], binnings['max(neEmEF)'] ])

    ax.set_xlabel(r'$\Delta \Phi_{jj}$')
    ax.set_ylabel(r'max(neEmEF)')
    fig.colorbar(im)

    # Figure title
    proc_to_title = {
        'VBF'   : r'VBF $H(inv)$',
        'Znunu' : r'$Z(\nu \nu)$',
        'Wlnu'  : r'$W(\ell \nu)$',
        'data'  : 'Data (2017F)'
    }
    ax.set_title(proc_to_title[proc])

    # Save figure
    outdir = f'./output/{tree_version}/2d_histos'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'{proc}_2d_histo.pdf')
    fig.savefig(outpath)
    print(f'MSG% File saved: {outpath}')

def main():
    # The ROOT file to use as an input for the 2D histogram
    input_root_dir = '/afs/cern.ch/work/a/aakpinar/public/forZeynep/VBF_trees/09Jul20'
    proc = sys.argv[1]

    root_files_for_procs = {
        'VBF'   : pjoin(input_root_dir, 'tree_VBF_HToInvisible_M125_pow_pythia8_2017.root'),
        'Znunu' : pjoin(input_root_dir, 'tree_ZJetsToNuNu_HT-400To600-mg_new_pmx_2017.root'),
        'Wlnu'  : pjoin(input_root_dir, 'tree_WJetsToLNu_HT-400To600-MLM_2017.root'),
        'data'  : pjoin(input_root_dir, 'tree_MET_2017F.root')
    }
    if proc != 'all':
        processes = [proc]
    else:
        # Plot for all processes 
        processes = root_files_for_procs.keys()

    # Get the tree version
    tree_version = os.path.basename(input_root_dir)

    for proc in processes:
        input_root_file = root_files_for_procs[proc]
        plot_2d_histogram(input_root_file, proc, tree_version)

if __name__ == '__main__':
    main()
