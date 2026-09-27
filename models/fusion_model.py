import torch.nn as nn
from .chemical_encoder import ChemicalEncoder
from .protein_encoder import ProteinEncoder
from .geometry_encoder import GeometryEncoder
from .interaction_encoder import InteractionEncoder
class ProtacGeoModel(nn.Module):
    def __init__(self, chemical_dim, protein_dim, geometry_dim, hidden=256, dropout=.2):
        super().__init__(); self.chem=ChemicalEncoder(chemical_dim,hidden); self.prot=ProteinEncoder(protein_dim,hidden); self.geo=GeometryEncoder(geometry_dim,128); self.inter=InteractionEncoder(hidden)
        self.head=nn.Sequential(nn.Linear(hidden+128,hidden),nn.LayerNorm(hidden),nn.SiLU(),nn.Dropout(dropout),nn.Linear(hidden,2))
    def forward(self, chemical,target,e3,geometry):
        c=self.chem(chemical); t=self.prot(target); e=self.prot(e3); g=self.geo(geometry); z=self.inter(t,e,c); return self.head(__import__('torch').cat([z,g],dim=-1))
