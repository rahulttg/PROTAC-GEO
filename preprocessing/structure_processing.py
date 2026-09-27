"""Real PDB parsing and geometry calculations for PROTAC-GEO."""
from __future__ import annotations
import math
from collections import defaultdict
from typing import Dict, List, Tuple
import numpy as np

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE"}
WATER = {"HOH","WAT","DOD"}

def _dist(a,b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b)))

def parse_pdb_text(text: str):
    atoms=[]
    for line in text.splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        if len(line) < 54: continue
        try:
            rec=line[0:6].strip(); atom=line[12:16].strip(); res=line[17:20].strip(); chain=line[21].strip() or "_"
            resseq=int(line[22:26].strip()); icode=line[26].strip()
            x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
        except ValueError:
            continue
        if res in WATER: continue
        elem=(line[76:78].strip() if len(line)>=78 else "") or atom[0]
        atoms.append({"record":rec,"atom":atom,"resname":res,"chain":chain,"resseq":resseq,"icode":icode,"x":x,"y":y,"z":z,"element":elem.upper()})
    return atoms

def _centroid(atoms):
    if not atoms: return None
    return tuple(np.mean(np.array([[a['x'],a['y'],a['z']] for a in atoms]),axis=0).tolist())

def analyze_pdb(text: str, radius=8.0):
    atoms=parse_pdb_text(text)
    if not atoms: raise ValueError("No ATOM/HETATM coordinates were found in the PDB text.")
    chains=defaultdict(list); residues=defaultdict(list); lig=[]
    for a in atoms:
        if a["record"]=="ATOM": chains[a["chain"]].append(a)
        key=(a["chain"],a["resseq"],a["icode"],a["resname"])
        if a["record"]=="ATOM": residues[key].append(a)
        elif a["resname"] not in AA3: lig.append(a)
    chain_names=sorted(chains)
    chain_centroids={c:_centroid(v) for c,v in chains.items()}
    ca_by_chain={c:[a for a in v if a['atom']=="CA"] for c,v in chains.items()}
    lys=[]
    for key,ras in residues.items():
        if key[3] != "LYS": continue
        ca=next((a for a in ras if a['atom']=="CA"), None)
        if ca is None: continue
        c=key[0]
        # A geometry-based exposure proxy: count nearby CA atoms in the same chain.
        coords=np.array([[a['x'],a['y'],a['z']] for a in ca_by_chain[c]],dtype=float)
        p=np.array([ca['x'],ca['y'],ca['z']])
        if len(coords):
            d=np.sqrt(((coords-p)**2).sum(axis=1)); neighbor_count=int(((d>1e-6)&(d<=10.0)).sum())
        else: neighbor_count=0
        lys.append({"chain":c,"residue":f"K{key[1]}{key[2]}","resseq":key[1],"x":ca['x'],"y":ca['y'],"z":ca['z'],"neighbor_count":neighbor_count})
    # Exposed = below-median local packing count. This is explicitly a proxy, not SASA.
    if lys:
        med=float(np.median([x['neighbor_count'] for x in lys]))
        for x in lys: x['exposure_proxy']=max(0.0,min(1.0,1.0-(x['neighbor_count']/(med+1.0))))
        exposed=[x for x in lys if x['neighbor_count']<=med]
    else:
        med=0.0; exposed=[]
    ligand_centroid=_centroid(lig)
    nearest=[]
    if ligand_centroid:
        for x in lys:
            x['distance_to_ligand_A']=_dist((x['x'],x['y'],x['z']),ligand_centroid)
        nearest=sorted(lys,key=lambda x:x['distance_to_ligand_A'])
        accessible=[x for x in lys if x.get('distance_to_ligand_A',1e9)<=radius and x['neighbor_count']<=med]
    else:
        accessible=exposed
    interface=None
    if len(chain_names)>=2:
        c1,c2=chain_names[:2]
        a1=np.array([[a['x'],a['y'],a['z']] for a in chains[c1]])
        a2=np.array([[a['x'],a['y'],a['z']] for a in chains[c2]])
        if len(a1) and len(a2):
            # blockwise min distances to avoid a huge O(N^2) temporary array
            contacts=0; min_d=999.0
            for start in range(0,len(a1),500):
                block=a1[start:start+500]
                dd=np.sqrt(((block[:,None,:]-a2[None,:,:])**2).sum(axis=2))
                contacts += int((dd<=5.0).sum()); min_d=min(min_d,float(dd.min()))
            cc1=chain_centroids[c1]; cc2=chain_centroids[c2]
            interface={"chains":[c1,c2],"contact_atoms_5A":contacts,"minimum_atom_distance_A":min_d,"centroid_distance_A":_dist(cc1,cc2)}
    return {
        "atom_count":len(atoms),"protein_atom_count":sum(1 for a in atoms if a['record']=="ATOM"),
        "hetero_atom_count":sum(1 for a in atoms if a['record']=="HETATM"),"chains":chain_names,
        "chain_centroids":chain_centroids,"ligand_centroid":ligand_centroid,
        "lysines":lys,"lysine_count":len(lys),"exposed_lysine_count":len(exposed),
        "accessible_lysine_count":len(accessible),"accessible_lysines":sorted(accessible,key=lambda x:x.get('distance_to_ligand_A',999.0))[:20],
        "nearest_lysines":nearest[:10],"interface":interface,
        "radius_A":radius,"exposure_method":"10 A CA-neighborhood packing proxy"
    }
