"""Minimal PyTorch training loop for the PROTAC-GEO multimodal model."""
import torch

def train_one_epoch(model, loader, optimizer, criterion, device="cpu"):
    model.train(); total=0.0
    for batch in loader:
        batch={k:v.to(device) for k,v in batch.items()}; y=batch.pop("labels"); optimizer.zero_grad(); logits=model(**batch); loss=criterion(logits,y); loss.backward(); optimizer.step(); total += loss.item()*len(y)
    return total/len(loader.dataset)
