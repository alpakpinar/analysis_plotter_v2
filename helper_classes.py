#!/usr/bin/env python

class Style:
    def __init__(self):
        # List of x-labels for each variable
        self.xlabels = {
            'mjj' : r'$M_{jj} \ (GeV)$',
            'dPhiLeadingJetMet'  : r'$\Delta \Phi (jet0, MET)$',
            'dPhiTrailingJetMet' : r'$\Delta \Phi (jet1, MET)$',
            'dPhiMoreCentralJetMet' : r'$\Delta \Phi (central jet, MET)$',
            'dPhiMoreForwardJetMet' : r'$\Delta \Phi (forward jet, MET)$',
            'mindPhiJetMet' : r'$min\Delta \Phi (jet, MET)$',
            'maxdPhiJetMet' : r'$max\Delta \Phi (jet, MET)$',
            'leadak4_eta'   : r'Leading Jet $\eta$',
            'trailak4_eta'  : r'Trailing Jet $\eta$',
            'leadak4_neHEF'   : 'Leading Jet Neutral Hadron Fraction',
            'trailak4_neHEF'  : 'Trailing Jet Neutral Hadron Fraction',
            'leadak4_neEmEF'  : 'Leading Jet Neutral EM Fraction',
            'trailak4_neEmEF' : 'Trailing Jet Neutral EM Fraction'
        }

        # List of figure titles for each region
        self.fig_titles = {
            'region A' : r'$\Delta \Phi_{jj} > 1.5$ & $max(neEmEF) < 0.8$',
            'region B' : r'$\Delta \Phi_{jj} > 1.5$ & $max(neEmEF) > 0.8$',
            'region C' : r'$\Delta \Phi_{jj} < 1.5$ & $max(neEmEF) > 0.8$',
            'region D' : r'$\Delta \Phi_{jj} < 1.5$ & $max(neEmEF) < 0.8$',
            'signal'   : r'$\Delta \Phi_{jj} < 1.5$ (SR Selection)',
            'dphijj_largerThan_1_5' : r'$\Delta \Phi_{jj} > 1.5$',
            'noCuts' : 'No Additional Cuts'
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
    def __init__(self):
        '''
        Create and store a dictionary mapping the regions to the cuts that are being used for the region.
        While calling this class, one should provide two variables, two low limits and high limits for each variable
        for cutting. One of the limits can be None.
        '''
        # For the ABCD method, there should be two variables being used
        # assert len(variables)   == 2
        # assert len(low_limits)  == 2
        # assert len(high_limits) == 2
        # self.variables   = variables
        # self.low_limits  = low_limits
        # self.high_limits = high_limits

        self.selections_by_region = {
            # Regions for ABCD method
            'region A' : [
                {'variable' : 'dphijj', 'low' : 1.5, 'high' : None},
                {'variable' : 'max(neEmEF)', 'low' : None, 'high' : 0.8}
            ],
            'region B' : [
                {'variable' : 'dphijj', 'low' : 1.5, 'high' : None},
                {'variable' : 'max(neEmEF)', 'low' : 0.8, 'high' : None}
            ],
            'region C' : [
                {'variable' : 'dphijj', 'low' : None, 'high' : 1.5},
                {'variable' : 'max(neEmEF)', 'low' : 0.8, 'high' : None}
            ],
            'region D' : [
                {'variable' : 'dphijj', 'low' : None, 'high' : 1.5},
                {'variable' : 'max(neEmEF)', 'low' : None, 'high' : 0.8}
            ],
            # Signal region selection: dphijj < 1.5
            'signal': [
                {'variable' : 'dphijj', 'low' : None, 'high' : 1.5}
            ],
            # Orthogonal to SR, dphijj > 1.5
            'dphijj_largerThan_1_5' : [
                {'variable' : 'dphijj', 'low' : 1.5, 'high' : None}
            ]
        }
