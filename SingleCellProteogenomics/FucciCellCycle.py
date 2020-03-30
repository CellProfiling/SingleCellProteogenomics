# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:26:52 2020

@author: antho
"""

class FucciCellCycle:
    '''
    Object representing the length of the FUCCI cell cycle phase transitions, which
    were manually determined by Diana M.
    '''
    def __init__(self):
        # Length of the cell cycle observed for the FUCCI cell line
        self.G1_LEN = 10.833 #hours (plus 10.833, so 13.458hrs for the S/G2 cutoff)
        self.G1_S_TRANS = 2.625 #hours (plus 10.833 and 2.625 so 25.433 hrs for the G2/M cutoff)
        self.S_G2_LEN = 11.975 #hours (this should be from the G2/M cutoff above to the end)
        self.M_LEN = 0.5 # We excluded M-phase from this analysis

        self.TOT_LEN = self.G1_LEN+self.G1_S_TRANS+self.S_G2_LEN

        self.G1_PROP = self.G1_LEN / self.TOT_LEN
        self.G1_S_PROP = self.G1_S_TRANS / self.TOT_LEN + self.G1_PROP
        self.S_G2_PROP = self.S_G2_LEN / self.TOT_LEN + self.G1_S_PROP