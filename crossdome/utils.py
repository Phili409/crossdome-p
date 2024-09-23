# Import needed libraries for this program 
import pandas as pd 
import os 
from typing import List 

def load_hla_database(file_path:str) -> pd.DataFrame:
    """
    Loads the HLA database (peptide dataset) from a CSV file into a pandas DataFrame.
    
    :param file_path: The file path of the CSV file containing the HLA database.
    :return: A pandas DataFrame with the HLA peptide data.
    """
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} does not exist.")
    
    try:
        # Attempt to open & load file 
        df:pd.DataFrame = pd.read_csv(file_path)
        # Return the dataframe 
        return df
    except Exception as e:
        # Error while attempting to open & load file 
        raise IOError(f"Error loading file {file_path}: {e}")
    
    

def load_background_peptides(df:pd.DataFrame, peptide_column:str = 'peptide') -> List[str]:
    """
    Extracts the background peptides from a pandas DataFrame.
    
    :param df: The pandas DataFrame containing the peptide data.
    :param peptide_column: The column name in the DataFrame containing the peptides.
    :return: A list of peptides from the DataFrame.
    """
    
    if peptide_column not in df.columns:
        raise KeyError(f"Column {peptide_column} not found in the DataFrame.")
    
    peptides: List[str] = df[peptide_column].tolist()
    return peptides



def save_results_to_csv(result_df:pd.DataFrame, file_path:str) -> None:
    """
    Saves a pandas DataFrame to a CSV file.
    
    :param result_df: The pandas DataFrame containing the results.
    :param file_path: The file path where the CSV should be saved.
    :return: None
    """
    
    try:
        result_df.to_csv(file_path, index=False)
        print(f"Results successfully saved to {file_path}")
    except Exception as e:
        raise IOError(f"Error saving file {file_path}: {e}")



def validate_peptide_length(peptides:List[str], expected_length:int = 9) -> List[str]:
    """
    Validates that all peptides in the list are of the expected length.
    
    :param peptides: A list of peptide strings.
    :param expected_length: The expected length of the peptides (default is 9-mers).
    :return: A list of valid peptides.
    :raises ValueError: If any peptide is not of the expected length.
    """
    
    invalid_peptides = [p for p in peptides if len(p) != expected_length]
    
    if invalid_peptides:
        raise ValueError(f"Some peptides are not of length {expected_length}: {invalid_peptides}")
    
    return peptides