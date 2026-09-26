import torch

def predict(model, batch):
    model.eval()
    with torch.no_grad():
        logits=model(**batch); probs=torch.softmax(logits,dim=-1); return probs[:,1]
