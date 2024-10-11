# Import needed libraries & packages 
import unittest 
import numpy as np 
import pandas as pd 

# Import core objects
from crossdome.core_classes import xrBackground, xrResult

# Import core functions 
from crossdome.core_functions import cross_compose, calculate_relatedness, cross_pair_summary, cross_write, cross_substitution_matrix

"""
Unit tests for the core functions in the CrossDome project.
These tests cover functions related to peptide comparisons and relatedness score calculations.

Tested Functions:
    - cross_compose: Compares a query peptide against a background set of peptides.
    - cross_pair_summary: Provides a summary of relatedness scores for a query and background peptides.
    - calculate_relatedness: Calculates the relatedness score between two peptides.

The tests ensure that these functions behave as expected and return valid results.
"""

class TestCoreFunctions (unittest.TestCase):
    """
    Unit tests for the core functions in CrossDome.
    
    Methods:
    
    setUp(): Prepares data for use in the tests.
    test_cross_compose(): Tests the cross_compose function to ensure proper output.
    test_cross_pair_summary(): Tests the cross_pair_summary function for correctness.
    test_calculate_relatedness(): Tests the calculate_relatedness function for valid scores.
    """
    
    def setUp (self):
        """
        Setup function to prepare data before each test runs.
        Initializes a query peptide and a list of background peptides for testing.
        """
        
        self.query:str           = "EVDPIGHLY"
        self.background_peptides = ["ESDPIVAQY", "EVDPIGHFY", "EVDPIGLLY"]
        self.background          = xrBackground(allele="HLA-A*01:01", peptides=self.background_peptides)
    
    def test_cross_compose (self):
        """
        Tests the cross_compose function to ensure that:
        
        - It returns an xrResult object with the correct query value.
        - The number of results matches the number of background peptides.
        - The 'relatedness_score' column is present in the result DataFrame.
        """
        
        result = cross_compose(self.query, self.background)
        self.assertEqual(result.query, self.query)
        self.assertEqual(len(result.result), len(self.background_peptides))
        self.assertTrue("relatedness_score" in result.result.columns)
        
    def test_cross_pair_summary(self):
        """
        Tests the cross_pair_summary function to verify:
        - The returned xrResult object has the correct query.
        - The number of rows matches the number of background peptides.
        - The 'relatedness_score' column is included in the result DataFrame.
        """
        
        result = cross_pair_summary (self.query, self.background)
        self.assertEqual(result.query, self.query)
        self.assertEqual(len(result.result), len(self.background_peptides))
        self.assertTrue("relatedness_score" in result.result.columns)

    def test_calculate_relatedness (self):
        """
        Tests the calculate_relatedness function to ensure that:
        - The function returns a numeric score.
        - The score is within the expected range of 0.0 to 1.0.
        """
        
        score = calculate_relatedness(self.query, self.background_peptides[0])
        self.assertIsInstance(score, float)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

if __name__ == '__main__':
    unittest.main()