def test_imports():
    from features.molecular_features import featurize_smiles
    from preprocessing.structure_processing import parse_pdb_atoms
    from evaluation.metrics import classification_metrics
