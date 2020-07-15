#!/usr/bin/env python

import re

class Style:
    def __init__(self):
        # List of x-labels for each variable
        self.xlabels = {
            'mjj' : r'$M_{jj} \ (GeV)$',
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
            'thirdJet_phi' : r'Third Jet $\phi$'
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

class Selection:
    def __init__(self, variables, thresholds, apply_recoil_cut=False, apply_jet_eta_cut=False, apply_jet_dphi_cut=False):
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

        first_variable, second_variable = self.variables
        first_thresh, second_thresh = self.thresholds

        self.selections_by_region = {
            # Regions for ABCD method
            'region A' : [
                {'variable' : first_variable, 'low' : first_thresh, 'high' : None},
                {'variable' : second_variable, 'low' : None, 'high' : second_thresh}
            ],
            'region B' : [
                {'variable' : first_variable, 'low' : first_thresh, 'high' : None},
                {'variable' : second_variable, 'low' : second_thresh, 'high' : None}
            ],
            'region C' : [
                {'variable' : first_variable, 'low' : None, 'high' : first_thresh},
                {'variable' : second_variable, 'low' : second_thresh, 'high' : None}
            ],
            'region D' : [
                {'variable' : first_variable, 'low' : None, 'high' : first_thresh},
                {'variable' : second_variable, 'low' : None, 'high' : second_thresh}
            ],
            # Signal region selection: dphijj < 1.5
            'signal': [
                {'variable' : first_variable, 'low' : None, 'high' : first_thresh}
                # {'variable' : 'dphijj', 'low' : None, 'high' : 1.5}
            ],
            # Orthogonal to SR, dphijj > 1.5
            'dphijj_largerThan_1_5' : [
                {'variable' : first_variable, 'low' : first_thresh, 'high' : None}
            ]
        }

        # Apply additional recoil>250 GeV cut if requested
        if apply_recoil_cut:
            for region in self.selections_by_region.keys():
                recoil_cut = {'variable' : 'recoil_pt', 'low' : 250, 'high' : None}
                self.selections_by_region[region].append(recoil_cut)

        if apply_jet_eta_cut:
            for region in self.selections_by_region.keys():
                jet_eta_cut = {'variable' : 'leadak4_trailak4_eta', 'low' : None, 'high' : 2.5}
                self.selections_by_region[region].append(jet_eta_cut)

        if apply_jet_dphi_cut:
            for region in self.selections_by_region.keys():
                jet_eta_cut = {'variable' : 'dPhi_TkMET_PFMET', 'low' : None, 'high' : 1.0}
                self.selections_by_region[region].append(jet_eta_cut)

        # Selection tag for output saving
        # Cleanup the dots/parantheses
        first_thresh_tag = re.sub('\.', '_', str(first_thresh))
        second_thresh_tag = re.sub('\.', '_', str(second_thresh))
        first_variable_tag = re.sub('\(|\)', '', first_variable)
        second_variable_tag = re.sub('\(|\)', '', second_variable)
        self.selection_tag = f'{first_variable_tag}_{first_thresh_tag}_{second_variable_tag}_{second_thresh_tag}'

        variable_to_fig_title = {
            'dphijj'      : r'$\Delta \Phi_{jj}$',
            'leadak4_neEmEF' : r'Lead AK4 neutral EM EF',
            'max(neEmEF)' : r'$max(neEmEF)$',
            'dPhi_TkMET_PFMET' : r'$\Delta \Phi (TkMET, PFMET)$',
        }

        # List of figure titles for each region
        if self.variables == ['dphijj', 'max(neEmEF)']:
            self.fig_titles = {
                'region A' : r'$\Delta \Phi_{{jj}} > {}$ & $max(neEmEF) < {}$'.format(first_thresh, second_thresh),
                'region B' : r'$\Delta \Phi_{{jj}} > {}$ & $max(neEmEF) > {}$'.format(first_thresh, second_thresh),
                'region C' : r'$\Delta \Phi_{{jj}} < {}$ & $max(neEmEF) > {}$'.format(first_thresh, second_thresh),
                'region D' : r'$\Delta \Phi_{{jj}} < {}$ & $max(neEmEF) < {}$'.format(first_thresh, second_thresh),
                'signal'   : r'$\Delta \Phi_{{jj}} < {}$ (SR Selection)'.format(first_thresh),
                'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
                'noCuts' : 'No Additional Cuts'
            }

        elif self.variables == ['dphijj', 'leadak4_neEmEF']:
            self.fig_titles = {
                'region A' : r'$\Delta \Phi_{{jj}} > {}$ & $leadAK4(neEmEF) < {}$'.format(first_thresh, second_thresh),
                'region B' : r'$\Delta \Phi_{{jj}} > {}$ & $leadAK4(neEmEF) > {}$'.format(first_thresh, second_thresh),
                'region C' : r'$\Delta \Phi_{{jj}} < {}$ & $leadAK4(neEmEF) > {}$'.format(first_thresh, second_thresh),
                'region D' : r'$\Delta \Phi_{{jj}} < {}$ & $leadAK4(neEmEF) < {}$'.format(first_thresh, second_thresh),
                'signal'   : r'$\Delta \Phi_{{jj}} < {}$ (SR Selection)'.format(first_thresh),
                'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
                'noCuts' : 'No Additional Cuts'
            }
            
        elif self.variables == ['dphijj', 'dPhi_TkMET_PFMET']:
            self.fig_titles = {
                'region A' : r'$\Delta \Phi_{{jj}} > {}$ & $\Delta \Phi (TkMET, PFMET) < {}$'.format(first_thresh, second_thresh),
                'region B' : r'$\Delta \Phi_{{jj}} > {}$ & $\Delta \Phi (TkMET, PFMET) > {}$'.format(first_thresh, second_thresh),
                'region C' : r'$\Delta \Phi_{{jj}} < {}$ & $\Delta \Phi (TkMET, PFMET) > {}$'.format(first_thresh, second_thresh),
                'region D' : r'$\Delta \Phi_{{jj}} < {}$ & $\Delta \Phi (TkMET, PFMET) < {}$'.format(first_thresh, second_thresh),
                'signal'   : r'$\Delta \Phi_{{jj}} < {}$ (SR Selection)'.format(first_thresh),
                'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
                'noCuts' : 'No Additional Cuts'
            }


