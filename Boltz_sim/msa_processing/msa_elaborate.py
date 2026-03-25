"""
Script to elaborate  msa and receive the correct numbering sequence and parameters for the pockets of BACE1, based on experimental results and analysis of the sequence used.

To be run after: 
> Having obtained the sequence from Uniprot
> Having run the notebook of Alphafold2 at https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb to generate msa
> Having run the scripts test_a3m.py and clean_msa.py to obtain a clean formatted msa.

What you should obtain from this script: 
> Relative numeration for the sequence
> Relative numeration for the pockets or the groups similar to the pockets 
"""