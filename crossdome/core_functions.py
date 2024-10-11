# Import needed libraries & packages
import numpy as np 
import pandas as pd
from typing import List, Dict
from scipy.stats import norm # for Z-scores & p-vals 
from datetime import datetime

# Import class objects 
from core_classes import xrBackground, xrResult

# Internal helper functions 
def _internal_checking_peptide (peptide:str) -> List[str]:
    """
    Checks if the peptide is a valid 9-mer and contains only standard amino acids 
    
    :param peptide: The peptide sequence 
    :return: Peptide split into a list of amino acids 
    :raises ValueError: If the peptide is not a valid 9-mer or contains non-standard amino acids.
    """
    
    if len(peptide) != 9:
        raise ValueError(f"CrossDome currently only supports 9-mer peptides, got size {len(peptide)} instead.")
    
    # Define a @set containing the standard amino acids 
    standard_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    
    if not all(residue in standard_amino_acids for residue in peptide):
        raise ValueError(f"Please check your input sequence {peptide}. Currently, only standard amino acids are supported.")

    return list(peptide) 
    
    
    
def _internal_related_distance (query_components:np.ndarray, subject_components:np.ndarray, position_weight:List[float]=None) -> float:
    """ 
    Calculates the relatedness score between two peptide components
    
    :param query_components: Matrix representation of query peptide 
    :param subject_components: Matrix representation of the subject peptide
    :param position_weight: Weights for each position in the peptide 
    :return: The relatedness score 
    """
    
    product_components = (query_components - subject_components) ** 2
    relatedness_score = np.sqrt(np.sum(product_components, axis=1))
    
    if position_weight is not None:
        relatedness_score = relatedness_score * np.sqrt(position_weight)
        
    relatedness_score = np.sum(relatedness_score) / query_components.shape[1]
    return relatedness_score



def _internal_percentile_rank (scores:List[float]) -> List[float]:
    """
    Calculates the percentile rank of each score.

    :param scores: A list of scores.
    :return: A list of percentile ranks for each score.
    """
    
    ranks:List = [(rank / (len(scores) - 1)) * 100 for rank in range(len(scores))]
    return ranks  



def _internal_matrix_correlation (query_components:np.ndarray, subject_components:np.ndarray) -> np.ndarray:
    """
    Computes the correlation between two peptide matrices.

    :param query_components: Matrix representation of the query peptide.
    :param subject_components: Matrix representation of the subject peptide.
    :return: The correlation matrix.
    """
    
    correlation_matrix = np.corrcoef(query_components, subject_components)
    return correlation_matrix



def cross_compose (query:str, background:xrBackground, position_weight:List[float] = None) -> xrResult:
    """
    This function compares a query peptide to a background set of peptides
    and returns an xrResult object containing relatedness scores.

    :param query: The query peptide (must be a 9-mer).
    :param background: An xrBackground object containing peptides to compare against.
    :param position_weight: A list of position weights (optional).
    :return: An xrResult object with comparison results.
    """

    # Validate the query peptide
    query_peptide = _internal_checking_peptide(query)

    # Use default position weights if none provided
    position_weight = position_weight or [1.0] * 9

    # Store results
    results: List[Dict[str, float]] = []

    # Iterate through each peptide in the background
    for element in background.peptides:
        # Validate the background peptide
        background_peptide = _internal_checking_peptide(element)

        # Calculate relatedness score
        score = _internal_related_distance(np.array(query_peptide), np.array(background_peptide), position_weight)

        # Count positive matches and mismatches
        num_positive = sum(q == c for q, c in zip(query, element))
        num_negative = 9 - num_positive

        # Add results
        results.append({
            'query': query,
            'subject': element,
            'relatedness_score': score,
            'num_positive': num_positive,
            'num_negative': num_negative
        })

    # Convert results to pandas DataFrame
    result_dataframe = pd.DataFrame(results)

    # Calculate Z-scores, p-values, percentile ranks, and ranks
    result_dataframe['zscore'] = (result_dataframe['relatedness_score'] - result_dataframe['relatedness_score'].mean()) / result_dataframe['relatedness_score'].std()
    result_dataframe['pvalue'] = norm.cdf(result_dataframe['zscore'])
    result_dataframe['percentile_rank'] = _internal_percentile_rank(result_dataframe['relatedness_score'].tolist())
    result_dataframe['rank'] = result_dataframe.index + 1

    # Return xrResult object
    return xrResult(query=query, result=result_dataframe, allele=background.allele, position_weight=position_weight)



def calculate_relatedness (query:str, candidate:str, position_weight:List[float] = None) -> float:
    """
    A function to calculate the relatedness score between a query peptide and a candidate peptide.

    :param query: The query peptide sequence.
    :param candidate: The candidate peptide sequence.
    :param position_weight: A list of weights applied to positions (optional).
    :return: A numeric relatedness score.
    """
    
    query_peptide = _internal_checking_peptide(query)
    candidate_peptide = _internal_checking_peptide(candidate)

    # Use default position weights if none provided
    position_weight = position_weight or [1.0] * 9

    # Calculate relatedness score
    return _internal_related_distance(np.array(query_peptide), np.array(candidate_peptide), position_weight)



def cross_pair_summary (query:str, background:xrBackground) -> xrResult:
    """
    Summarizes the relatedness scores between the query and background peptides.

    :param query: The query peptide sequence.
    :param background: An xrBackground object.
    :return: An xrResult object containing the summary statistics.
    """
    
    # Get the results from cross_compose
    xr_result = cross_compose(query=query, background=background)

    # The result DataFrame already contains the necessary information
    return xr_result



def cross_write (result:xrResult, file_path:str) -> None:
    """
        Exports the result DataFrame to a CSV file.

        :param result: An xrResult object containing peptide comparison results.
        :param file_path: The file path to save the CSV.
        :return: None
    """
    
    result.result.to_csv(file_path, index=False)