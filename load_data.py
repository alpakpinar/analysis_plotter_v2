#!/usr/bin/env python

import re
import csv
import os
import sys
from pprint import pprint
import uproot
import numpy as np

pjoin = os.path.join

# Physics process to dataset regex mapping
dataset_mapping = {
    'Diboson'  : re.compile('tree_(WW|WZ|ZZ)_.*'),
    'Top'      : re.compile('tree_(TTJets-amcatnloFXFX|ST_((s|t)-channel|tW)_(anti)?top.*inclusiveDecays.*).*'),
    'EWKW'     : re.compile('tree_EWKW(Minus|Plus).*'),
    'EWKZLL'   : re.compile('tree_EWKZ.*ZToLL.*'),
    'EWKZNuNu' : re.compile('tree_EWKZ.*ZToNuNu.*')
}

# Several different eta binnings for TF calculation in ABCD method
eta_binnings = {
    'very_fine' : np.linspace(-5,5,26),
    'fine' : list(range(-5,6)),   
    'coarse' : list(range(-5,6,2))
}

# Several different abs(eta) binnings for TF calculation in ABCD method
abs_eta_binnings = {
    'very_fine' : np.linspace(0,5,26),
    'fine' : np.arange(0,6),   
    'coarse_largeEta' : list(np.arange(0,3.6,0.2)) + [5] # One transfer factor for 3.4 < |eta| < 5.0
}

def event_isin_hfhf_category(events):
    '''Return a mask representing if the event is in HF-HF category (if both the leading jets are in HF).'''
    leadak4_abseta = np.abs(events['leadak4_eta'].array())
    trailak4_abseta = np.abs(events['trailak4_eta'].array())
    
    event_in_hfhf_category = (leadak4_abseta > 3.0) & (trailak4_abseta > 3.0)
    return event_in_hfhf_category

def get_data_from_csv(csv_file):
    '''Load data from the input CSV file: XS+sumw for each MC'''
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        # Construct the dataset --> xs/sumw map
        d = {row[0] : float(row[1])/float(row[2]) for row in reader if 'Dataset' not in row}
    return d

def load_data(inpath, process, csv_file, variable, cuts=None, eta_binning='very_fine', jes_variation=None, apply_cleaning_cuts=False, veto_hfhf=False):
    '''
    From the given input path, load the weighted and scaled histograms as a function of 
    the requested variable. Use the selections provided in the selection_dict variable. 
    '''
    binning = {
        'mjj' : list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500],
        'detajj' : np.linspace(0,10,51), 
        'dphijj' : np.linspace(0,3,31),
        'thirdJet_pt'  : list(range(30,370,20)),
        'thirdJet_eta' : np.linspace(-5,5,21),
        'thirdJet_phi' : np.linspace(-3.5,3.5,36),
        'nJet' : np.arange(0,10),
        'HT_jetsInHF' : np.linspace(10,510,26),
        'HTmiss_jetsInHF_pt' : np.linspace(10,510,26),
        'max(neEmEF)' : np.arange(0,1.05,0.05)
    }
    # Get XS + sumw scaling factors
    xs_sumw_scale = get_data_from_csv(csv_file)
    lumi = 41500

    # Get the relevant files for this process
    filelist = [f for f in os.listdir(inpath) if f.endswith('.root')]
    if process in dataset_mapping.keys():
        r = dataset_mapping[process]
    else:
        r = re.compile(f'^tree_{process}.*root')
    files = [pjoin(inpath, f) for f in list(filter(r.match, filelist))]
    # In the filelist, look for duplicate files:
    # If new_pmx version is available for a dataset, use that one instead and drop the regular file
    files_to_use = []
    for f in files:
        # Use the new_pmx files or extensions
        if re.match('.*new_+pmx.*', f) or re.match('.*ext.*', f):
            files_to_use.append(f)
            continue
        tag = f.split('/')[-1].replace('_2017.root' , '')
        new_pmx_exists = pjoin(inpath, f'{tag}_new_pmx_2017.root') in files
        if new_pmx_exists:
            # Do not append this file to the filelist 
            continue
        
        files_to_use.append(f)

    # Combined histogram will be stored in this variable 
    histo = None

    # Regions to look for different JES variations, default is the central one
    regions_to_look = {
        'central' : 'sr_vbf',
        'jesUp' : 'sr_vbf_jesTotalUp',
        'jesDown' : 'sr_vbf_jesTotalDown',
    }

    # Loop over each file and read contents
    for filename in files_to_use:
        filetag = filename.split('/')[-1].replace('.root', '').replace('tree_', '')
        # print(f'MSG% Working on: {filetag}')
        try:
            if jes_variation is not None:
                if process != 'MET':
                    region_to_look = regions_to_look[jes_variation]
                else:
                    region_to_look = 'sr_vbf'
            else:
                region_to_look = 'sr_vbf'
            events = uproot.open(filename)[region_to_look]
        except KeyError:
            print(f'MSG% WARNING: Could not find events in: {filename}, skipping.')
            continue

        # Event selection
        if cuts is not None:
            # Initialize the mask with all True values
            if variable in ['absEta', 'max(neEmEF)']:
                mask = np.ones_like(events['leadak4_eta'].array(), dtype=bool)
            else:
                mask = np.ones_like(events[variable].array(), dtype=bool)
            # Update the mask with each selection content
            for cut in cuts:
                mask_update = cut.get_mask(events)
                mask = mask & mask_update
            # Apply the cleaning cuts as well (VecB, VecDPhi), if requested (by default it is not)
            # In 29Jul20 trees, these are readily saved as a boolean
            if apply_cleaning_cuts:
                cleaning_mask = events['pass_cleaning_cut'].array().astype(bool)
                mask = mask & cleaning_mask
            
            # Veto HF-HF events if requested
            if veto_hfhf:
                hfhf_cat = event_isin_hfhf_category(events)
                mask = mask & (~hfhf_cat)
                
        else:
            mask = np.ones_like(events[variable].array(), dtype=bool)

        # Extract the variable, store as a histogram and scale
        if variable == 'absEta':
            var = np.abs(events['leadak4_eta'].array() )[mask]
        elif variable == 'max(neEmEF)':
            leadak4_neEmEF = np.abs(events['leadak4_neEmEF'].array())[mask]
            trailak4_neEmEF = np.abs(events['trailak4_neEmEF'].array())[mask]
            var = np.maximum(leadak4_neEmEF, trailak4_neEmEF)
        else:
            var = events[variable].array()[mask]
        # Get the weights for MC
        weights = np.ones_like(var)
        if process != 'MET':
            keys = [key.decode('utf-8') for key in events.keys()]
            weight_names = [key for key in keys if re.match('^weight_(?!(ele|mu|photon)).*', key)]
            # print(f'MSG% Weights being used: {" ".join(weight_names)}')        
            for weight in weight_names:
                weights *= events[weight].array()[mask]

        # Construct the histogram
        if variable in binning.keys():
            bins = binning[variable]
        elif re.match('.*dPhi.*', variable):
            bins = np.linspace(0,3.5,51)
        # Get the eta binning to be used (could be different for ABCD method calculations)
        elif re.match('.*ak4.*eta.*', variable):
            bins = eta_binnings[eta_binning]
        elif variable == 'absEta':
            bins = abs_eta_binnings[eta_binning]
        # Binning for jet energy fractions
        elif re.match('.*ak4.*EF', variable):
            bins = np.linspace(0,1,51)
        # Number of primary vertices
        elif re.match('PV_npvs.*', variable):
            bins = np.arange(0,40)
        else:
            raise RuntimeError(f'No binning found for variable: {variable}')
        h, bins = np.histogram(var, bins=bins, weights=weights)
        
        # Scale the MC histogram (do not scale data)
        if process != 'MET':
            scale_factor = xs_sumw_scale[filetag] * lumi
            h = h * scale_factor

        if histo is not None:
            histo += h
        else:
            histo = h

    # Return the histogram and the bins
    return histo, bins
