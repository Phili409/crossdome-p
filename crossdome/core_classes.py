# Imported libraries
import pandas as pd
import datetime
from typing import List 



class XRBackground (object):
    """
        A class to hold background peptide data and statistics.

        Attributes:
        allele (str): The HLA allele associated with the peptides.
        peptides (List[str]): A list of 9-mer peptide sequences.
        stats (dict): A dictionary containing statistical information.
    """

    def __init__(self, allele:str, peptides:List[str], stats:dict):
        self.allele = allele
        self.peptides = peptides
        self.stats = stats
        self.validate_peptides()

    def validate_peptides(self):
        """
            Validates that all peptides are 9-mers.

            Raises:
            ValueError: If any is not a 9-mer.
        """
        
        if not all(len(peptide) == 9 for peptide in self.peptides):
            raise ValueError("CrossDome only supports 9-mer peptides.")



class XRResult (object):
    """
        The main CrossDome result object.

        Attributes:
        query (str): The query peptide sequence.
        result (pd.DataFrame): DataFrame containing the relatedness ranking.
        allele (str): The HLA allele used in the analysis.
        expression (dict): Expression data.
        analysis (dict): Sequence and immunogenicity analysis data.
        position_weight (List[float]): Weight vector for peptide positions.
        timestamp (str): Timestamp of when the analysis was run.
    """

    def __init__(
        self,
        query: str,
        result: pd.DataFrame,
        allele: str,
        expression: dict = None,
        analysis: dict = None,
        position_weight: List[float] = None,
        timestamp: str = None,
    ):
        self.query = query
        self.result = result
        self.allele = allele
        self.expression = expression or {}
        self.analysis = analysis or {}
        self.position_weight = position_weight or [1.0] * 9
        self.timestamp = timestamp or datetime.now().strftime("%Y-%m-%d %H:%M:%S")
