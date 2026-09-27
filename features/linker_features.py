"""Linker descriptors. If an isolated linker SMILES is supplied, compute real descriptors."""
from .molecular_features import featurize_smiles

def linker_descriptor(linker_smiles: str):
    return featurize_smiles(linker_smiles, n_bits=64)
