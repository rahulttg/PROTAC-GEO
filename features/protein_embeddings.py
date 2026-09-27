"""ESM-2 embedding utilities. No random vectors are generated."""
import numpy as np

def esm2_embedding(sequence: str, model_name: str = "facebook/esm2_t6_8M_UR50D") -> np.ndarray:
    try:
        import torch
        from transformers import AutoTokenizer, AutoModel
    except ImportError as exc:
        raise RuntimeError("Install transformers and torch for ESM-2 embeddings.") from exc
    sequence = "".join(sequence.split()).upper()
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()
    with torch.no_grad():
        tokens = tokenizer(sequence, return_tensors="pt", truncation=True, max_length=1022)
        outputs = model(**tokens)
        mask = tokens["attention_mask"].unsqueeze(-1)
        hidden = outputs.last_hidden_state * mask
        pooled = hidden.sum(dim=1) / mask.sum(dim=1).clamp(min=1)
    return pooled.squeeze(0).cpu().numpy().astype(np.float32)
