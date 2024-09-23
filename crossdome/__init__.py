# Import modules and functions exposing them at package level
from .core_functions import cross_compose, cross_pair_summary, cross_substitution_matrix, cross_write
from .quant import peptide_similarity, overall_similarity, mismatch_distribution, peptide_distance
from .utils import load_hla_database, load_background_peptides, save_results_to_csv, validate_peptide_length
from .visualization import plot_similarity_heatmap, plot_mismatch_distribution, plot_relatedness_score_distribution, plot_peptide_expression

# Define the version of your package
__version__ = "1.0.0"

"""
CrossDome: A package for peptide comparison and analysis

This package includes modules for:
- Peptide similarity analysis
- Quantitative peptide comparison functions
- Visualization tools for peptide data
- Utility functions for loading and saving peptide datasets

Modules:
- core_functions: Core logic for peptide analysis
- quant: Quantitative analysis functions
- utils: Utility functions for file I/O and validation
- visualization: Functions for visualizing peptide-related data
"""