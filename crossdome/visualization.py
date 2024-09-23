# Import needed libraries 
import matplotlib.pyplot as plt 
import seaborn as sns 
import pandas as pd 
from typing import List 



def plot_similarity_heatmap(similarity_matrix:pd.DataFrame, title:str = "Peptide Similarity Heatmap") -> None:
    """
    Plots a heatmap of the similarity scores between peptides.
    
    :param similarity_matrix: A pandas DataFrame where each entry is a similarity score between peptides.
    :param title: The title of the heatmap plot.
    :return: None
    """
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(similarity_matrix, annot=True, cmap="coolwarm", linewidths=.5)
    plt.title(title)
    plt.show()



def plot_mismatch_distribution(mismatch_counts:List[int], title:str = "Mismatch Distribution") -> None:
    """
    Plots a bar plot showing the distribution of mismatch counts between peptides.
    
    :param mismatch_counts: A list of mismatch counts for each peptide comparison.
    :param title: The title of the bar plot.
    :return: None
    """
    
    plt.figure(figsize=(8, 6))
    plt.bar(range(len(mismatch_counts)), mismatch_counts, color='skyblue')
    plt.xlabel("Peptide Index")
    plt.ylabel("Number of Mismatches")
    plt.title(title)
    plt.show()



def plot_relatedness_score_distribution(scores:List[float], title:str = "Relatedness Score Distribution") -> None:
    """
    Plots a histogram of the relatedness scores between peptides.
    
    :param scores: A list of relatedness scores between peptides.
    :param title: The title of the histogram.
    :return: None
    """
    
    plt.figure(figsize=(8, 6))
    plt.hist(scores, bins=10, color='green', edgecolor='black', alpha=0.7)
    plt.xlabel("Relatedness Score")
    plt.ylabel("Frequency")
    plt.title(title)
    plt.show()



def plot_peptide_expression(expression_df:pd.DataFrame, peptide:str, title:str = "Peptide Expression Profile") -> None:
    """
    Plots the gene expression profile for a specific peptide across tissues.
    
    :param expression_df: A pandas DataFrame where each row represents a tissue and its expression level.
    :param peptide: The specific peptide whose expression profile is being plotted.
    :param title: The title of the bar plot.
    :return: None
    """
    
    if peptide not in expression_df.columns:
        raise KeyError(f"Peptide {peptide} not found in the DataFrame.")
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x=expression_df.index, y=expression_df[peptide], color='blue')
    plt.xlabel("Tissues")
    plt.ylabel("Expression Level")
    plt.title(title)
    plt.xticks(rotation=90)
    plt.show()