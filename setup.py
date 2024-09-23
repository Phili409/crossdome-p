# Import needed libraries 
from setuptools import setup, find_packages

# Define our setup 
setup(
    name="crossdome-p",  # Package Name
    version="1.0.0",     # Package Version
    author="Antunes Lab, Gabriel Galvez",  # Author name
    author_email="TBD",  # Author email
    description="CrossDome: Peptide comparison and analysis tool",
    long_description=open("README.md").read(),  # Detailed description from README.md
    long_description_content_type="text/markdown",
    url="https://github.com/Phili409/crossdome-p",  # GitHub URL of the project
    packages=find_packages(),  # Automatically find all packages
    install_requires=[  # Required dependencies
        "numpy",
        "pandas",
        "matplotlib",
        "seaborn"
    ],
    classifiers=[  # Categorization
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',  # Python version requirement
)
