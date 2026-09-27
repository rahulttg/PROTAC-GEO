import torch.nn as nn
class ChemicalEncoder(nn.Module):
    def __init__(self, in_dim, hidden=256):
        super().__init__(); self.net=nn.Sequential(nn.Linear(in_dim,hidden),nn.LayerNorm(hidden),nn.SiLU(),nn.Dropout(.2))
    def forward(self,x): return self.net(x)
