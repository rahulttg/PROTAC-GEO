"""Real RDKit molecular features for PROTAC-GEO."""
from __future__ import annotations
import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors, AllChem

def calculate_molecular_features(smiles: str, radius: int = 2, n_bits: int = 256):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Invalid SMILES string.")
    fp=AllChem.GetMorganFingerprintAsBitVect(mol,radius=radius,nBits=n_bits)
    arr=np.zeros(n_bits,dtype=np.int8); [arr.__setitem__(i,int(v)) for i,v in enumerate(fp)]
    f={
        "molecular_weight":float(Descriptors.MolWt(mol)),"logP":float(Crippen.MolLogP(mol)),"tpsa":float(rdMolDescriptors.CalcTPSA(mol)),
        "rotatable_bonds":int(Lipinski.NumRotatableBonds(mol)),"hbd":int(Lipinski.NumHDonors(mol)),"hba":int(Lipinski.NumHAcceptors(mol)),
        "rings":int(rdMolDescriptors.CalcNumRings(mol)),"aromatic_rings":int(Lipinski.NumAromaticRings(mol)),"fraction_csp3":float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "heavy_atoms":int(rdMolDescriptors.CalcNumHeavyAtoms(mol)),"formal_charge":int(Chem.GetFormalCharge(mol)),
        "exact_mass":float(Descriptors.ExactMolWt(mol)),"fingerprint_bits":int(arr.sum())
    }
    return f, arr.tolist()

def generate_conformer_metrics(smiles: str, n_confs: int = 12):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Invalid SMILES string.")
    mol=Chem.AddHs(mol)
    params=AllChem.ETKDGv3(); params.randomSeed=42
    ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=n_confs,params=params))
    if not ids: return {"conformer_count":0,"mmff_energies":[],"energy_spread_kcal_mol":None,"flexibility_proxy":0.0}
    energies=[]
    for cid in ids:
        try:
            ff=AllChem.MMFFGetMoleculeForceField(mol,AllChem.MMFFGetMoleculeProperties(mol),confId=cid)
            energies.append(float(ff.CalcEnergy()) if ff else float('nan'))
        except Exception: energies.append(float('nan'))
    finite=[e for e in energies if np.isfinite(e)]
    spread=float(max(finite)-min(finite)) if len(finite)>=2 else None
    # rotatable bond + conformer diversity proxy; not a thermodynamic entropy.
    rot=int(Lipinski.NumRotatableBonds(Chem.RemoveHs(mol)))
    flex=float(min(1.0,(rot/20.0)+(len(ids)/30.0)))
    return {"conformer_count":len(ids),"mmff_energies":finite[:12],"energy_spread_kcal_mol":spread,"flexibility_proxy":flex}
