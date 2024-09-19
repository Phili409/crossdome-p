# Import needed libraries
import numpy as np 
import pandas as pd
from typing import List, Dict

def cross_compose (query:str, background:List[str]) -> pd.DataFrame:
    """
        This function compares a query peptide to a background set of peptides
        and returns a dataframe of relatedness scores.
        
        :param query: The query peptide (must be a 9-mer).
        :param background: List of background peptides to compare against.
        :return: A pandas DataFrame with comparison results.
    """

    # Validate query, ensure a 9-mer query 
    if len(query) != 9:
        raise ValueError(f"Query peptide { query }, must be 9-mers got { len(query) } instead.")
    
    # Declare [results] variable, serves as datastructure to hold data 
    results:List[Dict[str, int]] = []
    
    # Iterate through each [element] _i in background peptides
    for element in background:
        # Calculate the relatedness score 
        score:float = calculate_relatedness(query, element)
        
        # Count positive mismatches & matches 
        num_positive:int = np.sum([q == c for (q, c) in zip(query, element)])
        num_negative:int = 9 - num_positive
        
        # Get parsed results into a dictionary 
        results.append({
            'query' : query,
            'subject' : element, 
            'relatedness score' : score,
            'num positive' : num_positive,
            'num negative' : num_negative
        })
        
    # Transform results to pandas DataFrame
    result_dataframe:pd.DataFrame = pd.DataFrame(results)
    return result_dataframe



def calculate_relatedness(query:str, candidate:str) -> float:
    """
    A simple function to calculate the relatedness score between a query peptide and a candidate peptide.
    
    :param query: The query peptide sequence.
    :param candidate: The candidate peptide sequence.
    :return: A numeric score representing relatedness.
    """
    
    return np.sum([q == c for q, c in zip(query, candidate)]) / len(query)



def cross_pair_summary(query:str, background:List[str]) -> pd.DataFrame:
    """
    Summarizes the relatedness scores between the query and background peptides.
    
    :param query: The query peptide sequence.
    :param background: A list of background peptides.
    :return: A pandas DataFrame with the relatedness score summary.
    """
    
    # Define results as a pandas dataframe
    result_dataframe:pd.DataFrame = cross_compose(query=query, background=background)
    
    # Statistics Summary: mean, standard deviation, minimum, maximum for relatedness scores 
    summary_stats:Dict[str, float] = {
        "mean relatedness" : result_dataframe["relatedness score"].mean(),
        "std relatedness" : result_dataframe["relatedness score"].std(),
        "min relatednes" : result_dataframe["relatedness score"].min(),
        "max relatednes" : result_dataframe["relatedness score"].max()
    }
    
    return pd.DataFrame([summary_stats])



def cross_substitution_matrix(query:str, candidate:str) -> List[int]:
    """
    Creates a matrix showing amino acid substitutions between the query and candidate peptides.
    
    :param query: The query peptide sequence.
    :param candidate: The candidate peptide sequence.
    :return: A list where 0 represents a match and 1 represents a mismatch.
    """
    
    # Define [substitution_matrix] via parsed list from query & candidate
    substitution_matrix:List[int] = [0 if q == c else 1 for q, c in zip(query, candidate)]
    return substitution_matrix

def cross_write(result_df:pd.DataFrame, file_path:str) -> None:
    """
    Exports the result DataFrame to a CSV file.
    
    :param result_df: The pandas DataFrame containing peptide comparison results.
    :param file_path: The file path to save the CSV.
    :return: None
    """
    
    result_df.to_csv(file_path, index=False)