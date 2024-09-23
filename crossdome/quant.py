# File library imports, for dependencies 
import numpy as np 
from typing import List 


def peptide_similarity(query:str, candidate:str) -> float:
    """
    Computes a similarity score between two peptides based on matching amino acids.
    
    :param query: The query peptide sequence (9-mer).
    :param candidate: The candidate peptide sequence (9-mer).
    :return: A float representing the similarity score (0 to 1).
    """
    
    if (len(query) != 9) or (len(candidate) != 9):
        # Base case 
        raise ValueError(f"Both query: {query} and candidate: {candidate} peptides must be 9-mers.")
    
    # Parse matches
    matches:int = np.sum([q == c for (q, c) in zip(query, candidate)])
    # Return matches
    return matches / len(query)



def peptide_mismatch_count(query:str, candidate:str) -> int:
    """
    Computes the number of mismatches between two peptide sequences.
    
    :param query: The query peptide sequence (9-mer).
    :param candidate: The candidate peptide sequence (9-mer).
    :return: An integer representing the number of mismatches.
    """
    
    if (len(query) != 9) or (len(candidate) != 9):
        # Base case 
        raise ValueError(f"Both query: {query} and candidate: {candidate} peptides must be 9-mers.")
    
    # Parse mismatches
    mismatches:int = np.sum([q != c for (q, c) in zip(query, candidate)])
    
    # Return mismatches
    return mismatches



def overall_similarity(query:str, background:List[str]) -> float:
    """
    Computes the overall average similarity score between a query peptide
    and a list of background peptides.
    
    :param query: The query peptide sequence (9-mer).
    :param background: A list of candidate peptide sequences.
    :return: A float representing the average similarity score.
    """
    
    # Parse similarities
    similarities:List[float] = [peptide_similarity(query, candidate) for candidate in background]
    # Calc avg for similarities
    average_similarity:float = np.mean(similarities)
    
    return average_similarity



def mismatch_distribution(query:str, background:List[str]) -> List[int]:
    """
    Computes the distribution of mismatches between a query peptide and a background of peptides.
    
    :param query: The query peptide sequence.
    :param background: A list of candidate peptide sequences.
    :return: A list of integers representing the mismatch counts for each peptide in the background.
    """
    
    mismatch_counts:List[int] = [peptide_mismatch_count(query, candidate) for candidate in background]
    return mismatch_counts



def peptide_distance(query:str, candidate:str) -> float:
    """
    Calculates a more advanced distance metric between two peptides, using a weighted score.
    
    :param query: The query peptide sequence.
    :param candidate: The candidate peptide sequence.
    :return: A float representing the weighted distance.
    """
    
    if len(query) != 9 or len(candidate) != 9:
        raise ValueError("Both query and candidate peptides must be 9-mers.")
    
    distances:List[float] = [1 if q != c else 0 for q, c in zip(query, candidate)]
    weighted_distance:float = np.sum(distances) 
    return weighted_distance