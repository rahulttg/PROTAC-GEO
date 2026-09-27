"""PROTAC SMILES validation and component-aware preparation."""
from typing import Dict

def canonicalize_smiles(smiles: str) -> str:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        return Chem.MolToSmiles(mol, canonical=True)
    except ImportError:
        return smiles.strip()

def validate_smiles(smiles: str) -> Dict[str, object]:
    result = {"valid": False, "canonical_smiles": None, "reason": None}
    try:
        result["canonical_smiles"] = canonicalize_smiles(smiles)
        result["valid"] = True
    except Exception as exc:
        result["reason"] = str(exc)
    return result
