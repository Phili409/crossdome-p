# Imported libraries
import pandas as pd
import datetime
from typing import List 



class xrBackground:
    """
    A class to represent the CrossDome background data holding peptides and related stats.
    
    Attributes:
        allele (str): The MHC Class I allele.
        peptides (List[str]): List of 9-mer peptides.
        stats (dict): A dictionary to hold statistics like off-target and database size.
    """
    
    def __init__(self, allele: str, peptides: List[str]):
        self.allele = allele
        self.peptides = peptides
        self.stats = {
            'off-target': 0,
            'database': len(peptides)
        }
        
        # Ensure all peptides are 9-mers
        if not all(len(pep) == 9 for pep in peptides):
            raise ValueError("All peptides must be 9-mers.")
    
    def __repr__(self):
        return f"xrBackground(allele={self.allele}, peptides_count={len(self.peptides)})"



class xrResult:
    """
    A class to represent the result of CrossDome peptide comparison.

    Attributes:
        query (str): The peptide target (9-mer).
        result (pd.DataFrame): The DataFrame containing relatedness ranking data.
        allele (str): The MHC Class I allele.
        expression (dict): A dictionary for storing expression data.
        analysis (dict): A dictionary for storing sequence and immunogenicity analysis.
        position_weight (List[float]): A list of numeric values derived from TCR hotspots.
        timestamp (str): The timestamp of the execution.
    """
    
    def __init__(self, query: str, result: pd.DataFrame, allele: str, position_weight: List[float]):
        self.query = query
        self.result = result
        self.allele = allele
        self.expression = {}
        self.analysis = {}
        self.position_weight = position_weight
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    def __repr__(self):
        return f"xrResult(query={self.query}, rank_count={len(self.result)})"
    
    def show(self):
        """Display the result DataFrame."""
        print(self.result)
    
    def add_expression(self, expression_data: dict):
        """Add expression data."""
        self.expression = expression_data
    
    def add_analysis(self, analysis_data: dict):
        """Add analysis data."""
        self.analysis = analysis_data