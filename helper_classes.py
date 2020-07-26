#!/usr/bin/env python

import re
import numpy as np

class Style:
    def __init__(self):
        # List of x-labels for each variable
        self.xlabels = {
            'mjj' : r'$M_{jj} \ (GeV)$',
            'detajj' : r'$\Delta \eta_{jj}$',
            'dphijj' : r'$\Delta \Phi_{jj}$',
            'dPhiLeadingJetMet'  : r'$\Delta \Phi (jet0, MET)$',
            'dPhiTrailingJetMet' : r'$\Delta \Phi (jet1, MET)$',
            'dPhiMoreCentralJetMet' : r'$\Delta \Phi (central jet, MET)$',
            'dPhiMoreForwardJetMet' : r'$\Delta \Phi (forward jet, MET)$',
            'dPhi_TkMET_PFMET' : r'$\Delta \Phi (TK MET, PF MET)$',
            'mindPhiJetMet' : r'$min\Delta \Phi (jet, MET)$',
            'maxdPhiJetMet' : r'$max\Delta \Phi (jet, MET)$',
            'leadak4_eta'   : r'Leading Jet $\eta$',
            'absEta'        : r'|$\eta$|',
            'trailak4_eta'  : r'Trailing Jet $\eta$',
            'leadak4_neHEF'   : 'Leading Jet Neutral Hadron Fraction',
            'trailak4_neHEF'  : 'Trailing Jet Neutral Hadron Fraction',
            'leadak4_neEmEF'  : 'Leading Jet Neutral EM Fraction',
            'trailak4_neEmEF' : 'Trailing Jet Neutral EM Fraction',
            'thirdJet_pt'  : r'Third Jet $p_T$',
            'thirdJet_eta' : r'Third Jet $\eta$',
            'thirdJet_phi' : r'Third Jet $\phi$',
            'PV_npvs'      : r'$N_{PV}$',
            'PV_npvsGood'  : r'Good $N_{PV}$',
            'nJet'  : r'$N_{jet}$',
            'HT_jetsInHF' : r'$H_T$ (Jets in HF)',
            'HTmiss_jetsInHF_pt' : r'$\vec{H}_{T,miss}$ (Jets in HF)',
            'max(neEmEF)' : 'Maximum neutral EM fraction'
        }

        # Pretty labels for legend for each process
        self.pretty_labels = {
            'ZJetsToNuNu' : r'QCD $Z\rightarrow \nu \nu$',
            'EWKZNuNu'    : r'EWK $Z\rightarrow \nu \nu$',
            'DYJetsToLL'  : r'QCD $Z\rightarrow \ell \ell$',
            'EWKZLL'      : r'EWK $Z\rightarrow \ell \ell$',
            'WJetsToLNu'  : r'QCD $W\rightarrow \ell \nu$',
            'EWKW'        : r'EWK $W\rightarrow \ell \nu$'
        }

class Cut:
    '''Helper class to represent a single cut, and retrieve the mask.'''
    def __init__(self, variable, low_thresh, high_thresh, special_apply=None, change_endcap_def=None):
        self.variable = variable
        self.low_thresh = low_thresh
        self.high_thresh = high_thresh
        self.special_apply = special_apply # Apply if one of the jets are in endcap? etc.
        self.change_endcap_def = change_endcap_def # Change maximum endcap coverage slightly, to |eta| = 3.2 from 3.0
        
    def get_mask(self, events):
        '''Create a boolean mask for the cut, reads the events object from outside.'''
        leadak4_abseta = np.abs(events['leadak4_eta'].array())
        trailak4_abseta = np.abs(events['trailak4_eta'].array())
        # Special case: Apply the cut only if one of the two leading jets is in endcap
        if self.special_apply == 'oneJetInEndcap':
            # If requested, slightly change the endcap maximum coverage 
            if self.change_endcap_def is not None:
                endcap_eta_largest = self.change_endcap_def
            else:
                endcap_eta_largest = 3.0 
            events_to_apply_cut = ((leadak4_abseta > 2.5) & (leadak4_abseta < endcap_eta_largest)) | ((trailak4_abseta > 2.5) & (trailak4_abseta < endcap_eta_largest))
        # Apply the cut only if none of the leading two jets is in HF
        elif self.special_apply == 'noJetInHF':
            events_to_apply_cut = (leadak4_abseta <= 3.0) & (trailak4_abseta <= 3.0)
        else:
            events_to_apply_cut = np.ones_like(leadak4_abseta, dtype=bool)

        # Handle specific variables that are not directly saved in the tree, e.g. max(neEmEF)
        if self.variable == 'max(neEmEF)':
            leadak4_neEmEF  = events['leadak4_neEmEF'].array()
            trailak4_neEmEF = events['trailak4_neEmEF'].array()
            variable_arr = np.maximum(leadak4_neEmEF, trailak4_neEmEF)
        elif self.variable in ['leadak4_trailak4_eta', 'more_central_leadingJet']:
            leadak4_eta  = np.abs(events['leadak4_eta'].array() )
            trailak4_eta = np.abs(events['trailak4_eta'].array() )
            variable_arr = np.minimum(leadak4_eta, trailak4_eta)
        elif self.variable == 'more_forward_leadingJet':
            leadak4_eta  = np.abs(events['leadak4_eta'].array() )
            trailak4_eta = np.abs(events['trailak4_eta'].array() )
            variable_arr = np.maximum(leadak4_eta, trailak4_eta)
        elif self.variable in ['absEta', 'leadak4_absEta']:
            variable_arr = np.abs(events['leadak4_eta'].array())
        elif self.variable == 'trailak4_absEta':
            variable_arr = np.abs(events['trailak4_eta'].array())
        else:
            variable_arr = events[self.variable].array()

        # Create the mask, depending on boundary conditions
        if (self.low_thresh is not None) and (self.high_thresh is not None):
            mask = (variable_arr > self.low_thresh) & (variable_arr < self.high_thresh)
        elif self.low_thresh is None:
            mask = variable_arr < self.high_thresh
        elif self.high_thresh is None:
            mask = variable_arr > self.low_thresh

        # Get the final mask: Apply the cut where wanted, for other events just return True
        final_mask = np.where(events_to_apply_cut, mask, True)
        return final_mask

class Selection:
    def __init__(self, variables, thresholds, apply_cuts=['recoil', 'met_dphi'], categorization=None):
        '''
        Create and store a dictionary mapping the regions to the cuts that are being used for the region.
        While calling this class, one should provide two variables, two low limits and high limits for each variable
        for cutting. One of the limits can be None.
        '''
        # For the ABCD method, there should be two variables being used
        assert len(variables)   == 2
        assert len(thresholds)  == 2
        self.variables   = variables
        self.thresholds  = thresholds

        self.first_variable, self.second_variable = self.variables
        self.first_thresh, self.second_thresh = self.thresholds

        self.apply_cuts = apply_cuts

        # Categorization for two leading jets, i.e. EE-HF. By default, no such categorization is applied.
        # If a specific categorization is chosen, additional two cuts on two leading jets will be applied.
        self.categorization = categorization

        self.get_selections_by_region()
        self.get_selection_tags()
        self.get_fig_titles()

    def get_cuts_for_category(self):
        '''If a categorization is specified, get the relevant cuts.'''
        # Jet eta cuts for each category
        categories_to_jet_eta_cut = {
            'Trk-Trk' : [
                Cut('more_central_leadingJet', low_thresh=None, high_thresh=2.5), 
                Cut('more_forward_leadingJet', low_thresh=None, high_thresh=2.5)
                ],
            'Trk-EE' : [
                Cut('more_central_leadingJet', low_thresh=None, high_thresh=2.5), 
                Cut('more_forward_leadingJet', low_thresh=2.5, high_thresh=3.0)
                ],
            'EE-EE' : [
                Cut('more_central_leadingJet', low_thresh=2.5, high_thresh=3.0), 
                Cut('more_forward_leadingJet', low_thresh=2.5, high_thresh=3.0)
                ],
            'Trk-HF' : [
                Cut('more_central_leadingJet', low_thresh=None, high_thresh=2.5), 
                Cut('more_forward_leadingJet', low_thresh=3.0, high_thresh=5.0)
                ],
            'EE-HF' : [
                Cut('more_central_leadingJet', low_thresh=2.5, high_thresh=3.0), 
                Cut('more_forward_leadingJet', low_thresh=3.0, high_thresh=5.0)
                ],
            'HF-HF' : [
                Cut('more_central_leadingJet', low_thresh=3.0, high_thresh=5.0), 
                Cut('more_forward_leadingJet', low_thresh=3.0, high_thresh=5.0)
                ],
        }
        return categories_to_jet_eta_cut[self.categorization]

    def get_selections_by_region(self):
        '''Get list of selections for each region.'''
        self.selections_by_region = {
            # Regions for ABCD method
            'region A' : [ 
                Cut(self.first_variable, low_thresh=self.first_thresh, high_thresh=None),
                Cut(self.second_variable, low_thresh=None, high_thresh=self.second_thresh)
            ],
            'region B' : [
                Cut(self.first_variable, low_thresh=self.first_thresh, high_thresh=None),
                Cut(self.second_variable, low_thresh=self.second_thresh, high_thresh=None)
            ],
            'region C' : [
                Cut(self.first_variable, low_thresh=None, high_thresh=self.first_thresh),
                Cut(self.second_variable, low_thresh=self.second_thresh, high_thresh=None)
            ],
            'region D' : [
                Cut(self.first_variable, low_thresh=None, high_thresh=self.first_thresh),
                Cut(self.second_variable, low_thresh=None, high_thresh=self.second_thresh)
            ],
            # Signal region selection: dphijj < 1.5
            'signal': [
                Cut(self.first_variable, low_thresh=None, high_thresh=self.first_thresh)
                # {'variable' : 'dphijj', 'low' : None, 'high' : 1.5}
            ],
            # Orthogonal to SR, dphijj > 1.5
            'dphijj_largerThan_1_5' : [
                Cut(self.first_variable, low_thresh=self.first_thresh, high_thresh=None)
            ]
        }

        # Additional cuts to be applied on all ABCD regions
        self.additional_cuts = {
            'recoil'   : Cut('recoil_pt', low_thresh=250, high_thresh=None),
            'max_neEmEF_0_85'   : Cut('max(neEmEF)', low_thresh=None, high_thresh=0.85),
            'dphijj_2_3'   : Cut('dphijj', low_thresh=None, high_thresh=2.3),
            'dphijj'   : Cut('dphijj', low_thresh=None, high_thresh=2.5),
            'absEta'   : Cut('absEta', low_thresh=3.0, high_thresh=None),
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
            'met_dphi_jetInEndcap_v2' : Cut('dPhi_TkMET_PFMET', low_thresh=None, high_thresh=0.75, special_apply='oneJetInEndcap', change_endcap_def=3.2),
            'leading_jet_pt_jetInEndcap_v2' : Cut('leadak4_pt', low_thresh=100, high_thresh=None, special_apply='oneJetInEndcap', change_endcap_def=3.2),
            'leading_jet_pt120_jetInEndcap_v2' : Cut('leadak4_pt', low_thresh=120, high_thresh=None, special_apply='oneJetInEndcap', change_endcap_def=3.2),
            # Cuts with endcap coverage slightly increased to |eta| 3.4
            'met_dphi_jetInEndcap_v3' : Cut('dPhi_TkMET_PFMET', low_thresh=None, high_thresh=0.75, special_apply='oneJetInEndcap', change_endcap_def=3.4),
            'leading_jet_pt_jetInEndcap_v3' : Cut('leadak4_pt', low_thresh=100, high_thresh=None, special_apply='oneJetInEndcap', change_endcap_def=3.4),
            'leading_jet_pt120_jetInEndcap_v3' : Cut('leadak4_pt', low_thresh=120, high_thresh=None, special_apply='oneJetInEndcap', change_endcap_def=3.4),
        }

        # Apply additional cuts if requested
        for cut_tag, cut in self.additional_cuts.items():
            if cut_tag in self.apply_cuts:
                for region in self.selections_by_region.keys():
                    self.selections_by_region[region].append(cut)

        if self.categorization is not None:
            categorization_cuts = self.get_cuts_for_category()
            for region in self.selections_by_region.keys():
                self.selections_by_region[region].extend(categorization_cuts)

    def get_selection_tags(self):
        '''
        Create selection tags based on the selections applied for ABCD method and the additional cuts applied.
        These tags will be used to name output directories to save plots made by applying the relevant selections.
        '''
        # Selection tag for output saving
        # Cleanup the dots/parantheses
        first_thresh_tag = re.sub('\.', '_', str(self.first_thresh))
        second_thresh_tag = re.sub('\.', '_', str(self.second_thresh))
        first_variable_tag = re.sub('\(|\)', '', self.first_variable)
        second_variable_tag = re.sub('\(|\)', '', self.second_variable)
        self.selection_tag = f'{first_variable_tag}_{first_thresh_tag}_{second_variable_tag}_{second_thresh_tag}'

        # Figure out which additional cuts are applied (to all ABCD regions)
        additional_cut_tags = []
        for cut_tag, cut in self.additional_cuts.items():
            if cut_tag in self.apply_cuts:
                additional_cut_tags.append(cut_tag)

        if len(additional_cut_tags) != 0:
            self.additional_selection_tag = '_'.join(additional_cut_tags)
        else:
            self.additional_selection_tag = 'noAdditionalSelection'

    def get_fig_titles(self):
        variable_to_fig_title = {
            'dphijj'      : r'$\Delta \Phi_{jj}$',
            'leadak4_neEmEF' : r'Lead AK4 neutral EM EF',
            'max(neEmEF)' : r'$max(neEmEF)$',
            'dPhi_TkMET_PFMET' : r'$\Delta \Phi (TkMET, PFMET)$',
        }

        # List of figure titles for each region
        if self.variables == ['dphijj', 'max(neEmEF)']:
            self.fig_titles = {
                'region A' : r'$\Delta \Phi_{{jj}} > {}$ & $max(neEmEF) < {}$'.format(self.first_thresh, self.second_thresh),
                'region B' : r'$\Delta \Phi_{{jj}} > {}$ & $max(neEmEF) > {}$'.format(self.first_thresh, self.second_thresh),
                'region C' : r'$\Delta \Phi_{{jj}} < {}$ & $max(neEmEF) > {}$'.format(self.first_thresh, self.second_thresh),
                'region D' : r'$\Delta \Phi_{{jj}} < {}$ & $max(neEmEF) < {}$'.format(self.first_thresh, self.second_thresh),
                'signal'   : r'$\Delta \Phi_{{jj}} < {}$ (SR Selection)'.format(self.first_thresh),
                'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
                'noCuts' : 'No Additional Cuts'
            }

        elif self.variables == ['dphijj', 'leadak4_neEmEF']:
            self.fig_titles = {
                'region A' : r'$\Delta \Phi_{{jj}} > {}$ & $leadAK4(neEmEF) < {}$'.format(self.first_thresh, self.second_thresh),
                'region B' : r'$\Delta \Phi_{{jj}} > {}$ & $leadAK4(neEmEF) > {}$'.format(self.first_thresh, self.second_thresh),
                'region C' : r'$\Delta \Phi_{{jj}} < {}$ & $leadAK4(neEmEF) > {}$'.format(self.first_thresh, self.second_thresh),
                'region D' : r'$\Delta \Phi_{{jj}} < {}$ & $leadAK4(neEmEF) < {}$'.format(self.first_thresh, self.second_thresh),
                'signal'   : r'$\Delta \Phi_{{jj}} < {}$ (SR Selection)'.format(self.first_thresh),
                'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
                'noCuts' : 'No Additional Cuts'
            }
            
        elif self.variables == ['dphijj', 'dPhi_TkMET_PFMET']:
            self.fig_titles = {
                'region A' : r'$\Delta \Phi_{{jj}} > {}$ & $\Delta \Phi (TkMET, PFMET) < {}$'.format(self.first_thresh, self.second_thresh),
                'region B' : r'$\Delta \Phi_{{jj}} > {}$ & $\Delta \Phi (TkMET, PFMET) > {}$'.format(self.first_thresh, self.second_thresh),
                'region C' : r'$\Delta \Phi_{{jj}} < {}$ & $\Delta \Phi (TkMET, PFMET) > {}$'.format(self.first_thresh, self.second_thresh),
                'region D' : r'$\Delta \Phi_{{jj}} < {}$ & $\Delta \Phi (TkMET, PFMET) < {}$'.format(self.first_thresh, self.second_thresh),
                'signal'   : r'$\Delta \Phi_{{jj}} < {}$ (SR Selection)'.format(self.first_thresh),
                'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
                'noCuts' : 'No Additional Cuts'
            }

        elif self.variables == ['dphijj', 'HT_jetsInHF']:
            self.fig_titles = {
                'region A' : r'$\Delta \Phi_{{jj}} > {}$ & $H_{{T,HF}} < {}$'.format(self.first_thresh, self.second_thresh),
                'region B' : r'$\Delta \Phi_{{jj}} > {}$ & $H_{{T,HF}} > {}$'.format(self.first_thresh, self.second_thresh),
                'region C' : r'$\Delta \Phi_{{jj}} < {}$ & $H_{{T,HF}} > {}$'.format(self.first_thresh, self.second_thresh),
                'region D' : r'$\Delta \Phi_{{jj}} < {}$ & $H_{{T,HF}} < {}$'.format(self.first_thresh, self.second_thresh),
                'signal'   : r'$\Delta \Phi_{{jj}} < {}$ (SR Selection)'.format(self.first_thresh),
                'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
                'noCuts' : 'No Additional Cuts'
            }

        elif self.variables == ['dphijj', 'HTmiss_jetsInHF_pt']:
            self.fig_titles = {
                'region A' : r'$\Delta \Phi_{{jj}} > {}$ & $\vec{{H}}_{{T,miss}} < {}$'.format(self.first_thresh, self.second_thresh),
                'region B' : r'$\Delta \Phi_{{jj}} > {}$ & $\vec{{H}}_{{T,miss}} > {}$'.format(self.first_thresh, self.second_thresh),
                'region C' : r'$\Delta \Phi_{{jj}} < {}$ & $\vec{{H}}_{{T,miss}} > {}$'.format(self.first_thresh, self.second_thresh),
                'region D' : r'$\Delta \Phi_{{jj}} < {}$ & $\vec{{H}}_{{T,miss}} < {}$'.format(self.first_thresh, self.second_thresh),
                'signal'   : r'$\Delta \Phi_{{jj}} < {}$ (SR Selection)'.format(self.first_thresh),
                'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
                'noCuts' : 'No Additional Cuts'
            }


