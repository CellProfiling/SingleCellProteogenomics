#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 17:13:54 2021

@author: anthony.cesnik
"""

import numpy as np
import unittest
from SingleCellProteogenomics import ProteinGaussianClustering, ProteinDataPreparation


class TestTest(unittest.TestCase):
    def setUp(self):
        protein_data = ProteinDataPreparation.ProteinData(False)

        # make example dataframes to zero center
        protein_data.green_fucci = np.array([1000, 2000, 3000, 2000])
        protein_data.red_fucci = np.array([1000, 2000, 1000, 2000])
        protein_data.u_plate = np.array(["1", "2"])
        protein_data.well_plate = np.array(["A5_1", "A6_1", "A5_2", "A5_2"])
        protein_data.plate = np.array(["1", "1", "2", "2"])

        self.protein = ProteinGaussianClustering.ProteinGaussianClustering(protein_data)
        self.protein.zero_center_fucci()

    def test_zero_centering(self):
        """Check that the log10 median is the center for each plate, for both red and green FUCCI colors"""
        for plate in self.protein.protein_data.u_plate:
            fromplate = self.protein.protein_data.plate == plate
            greenplatemedian = np.log10(
                np.median(self.protein.protein_data.green_fucci[fromplate])
            )
            redplatemedian = np.log10(
                np.median(self.protein.protein_data.red_fucci[fromplate])
            )
            assert all(
                np.log10(self.protein.protein_data.green_fucci[fromplate])
                == self.protein.log_green_fucci_zeroc[fromplate] + greenplatemedian
            )
            assert all(
                np.log10(self.protein.protein_data.red_fucci[fromplate])
                == self.protein.log_red_fucci_zeroc[fromplate] + redplatemedian
            )
            assert all(
                np.greater_equal(self.protein.log_green_fucci_zeroc_rescale, 0)
            )
            assert all(
                np.greater_equal(self.protein.log_red_fucci_zeroc_rescale, 0)
            )


if __name__ == "__main__":
    unittest.main()
