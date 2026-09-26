import torch
import torch.nn as nn
class InteractionEncoder(nn.Module):
    def __init__(self, hidden=256):
        super().__init__(); self.attn=nn.MultiheadAttention(hidden,4,batch_first=True); self.norm=nn.LayerNorm(hidden)
    def forward(self,target,e3,chemical):
        x=torch.stack([target,e3,chemical],dim=1); y,_=self.attn(x,x,x); return self.norm(y.mean(dim=1))
