# Import needed libraries & packages
import numpy as np 
import pandas as pd
from typing import List, Dict
from scipy.stats import norm # for Z-scores & p-vals 
from datetime import datetime

# Import class objects 
from core_classes import xrBackground, xrResult

# Import needed libraries
import numpy as np
import pandas as pd
from typing import List, Dict
from scipy.stats import norm
from core_classes import xrBackground, xrResult  # Import classes


def cross_compose(query:str, background:xrBackground, position_weight:List[float] = None) -> xrResult:
    """
        This function compares a query peptide to a background set of peptides
        and returns an xrResult object containing relatedness scores.

        :param query: The query peptide (must be a 9-mer).
        :param background: An xrBackground object containing peptides to compare against.
        :param position_weight: A list of position weights (optional).
        :return: An xrResult object with comparison results.
    """

    # Ensure that the query is a 9-mer
    if len(query) != 9:
        raise ValueError(f"Query peptide {query} must be 9-mers, got {len(query)} instead.")

    position_weight = position_weight or [1.0] * 9  # Default to equal weights if none provided

    # Store results
    results:List[Dict[str, float]] = []

    # Iterate through each peptide in the background
    for element in background.peptides:
        # Calculate relatedness score
        score:float = calculate_relatedness(query, element, position_weight)

        # Count positive matches & mismatches
        num_positive:int = np.sum([q == c for q, c in zip(query, element)])
        num_negative:int = 9 - num_positive

        # Add results
        results.append({
            'query': query,
            'subject': element,
            'relatedness_score': score,
            'num_positive': num_positive,
            'num_negative': num_negative
        })

    # Convert results to pandas DataFrame
    result_dataframe:pd.DataFrame = pd.DataFrame(results)
    
    result_dataframe['zscore']          = (result_dataframe['relatedness_score'] - result_dataframe['relatedness_score'].mean()) / result_dataframe['relatedness_score'].std()
    result_dataframe['pvalue']          = norm.cdf(result_dataframe['zscore'])
    result_dataframe['percentile_rank'] = result_dataframe['relatedness_score'].rank(pct=True)
    result_dataframe['rank']            = result_dataframe.index + 1

    # Return xrResult object
    return xrResult(query=query, result=result_dataframe, allele=background.allele, position_weight=position_weight)


def calculate_relatedness(query:str, candidate:str, position_weight:List[float] = None) -> float:
    """
        A simple function to calculate the relatedness score between a query peptide and a candidate peptide.

        :param query: The query peptide sequence.
        :param candidate: The candidate peptide sequence.
        :param position_weight: A list of weights applied to positions (optional).
        :return: A numeric relatedness score.
    """
    
    if position_weight is None:
        position_weight = [1.0] * 9  # Default weights

    # Calc relatedness score
    return np.sum([q == c for q, c in zip(query, candidate)]) / len(query)


def cross_pair_summary(query:str, background:xrBackground) -> xrResult:
    """
        Summarizes the relatedness scores between the query and background peptides.

        :param query: The query peptide sequence.
        :param background: An xrBackground object.
        :return: An xrResult object containing the summary statistics.
    """
    
    # Call cross_compose to parse comparison results
    result_dataframe:pd.DataFrame = cross_compose(query=query, background=background).result

    # Calc summary of stats
    summary_stats:Dict[str, float] = {
        "mean relatedness" : result_dataframe["relatedness_score"].mean(),
        "std relatedness"  : result_dataframe["relatedness_score"].std(),
        "min relatedness"  : result_dataframe["relatedness_score"].min(),
        "max relatedness"  : result_dataframe["relatedness_score"].max()
    }

    return xrResult(query=query, result=result_dataframe, allele=background.allele, position_weight=[1.0] * 9)


def cross_substitution_matrix(query:str, candidate:str) -> List[int]:
    """
        Creates a matrix showing amino acid substitutions between the query and candidate peptides.

        :param query: The query peptide sequence.
        :param candidate: The candidate peptide sequence.
        :return: A list where 0 represents a match and 1 represents a mismatch.
    """
    
    return [0 if q == c else 1 for q, c in zip(query, candidate)]


def cross_write(result:xrResult, file_path:str) -> None:
    """
        Exports the result DataFrame to a CSV file.

        :param result: An xrResult object containing peptide comparison results.
        :param file_path: The file path to save the CSV.
        :return: None
    """
    
    result.result.to_csv(file_path, index=False)