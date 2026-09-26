"""Optional ESM-2 embedding helper.
Install fair-esm separately and provide a local model checkpoint.
The core PROTAC-GEO web app does not download large model weights automatically.
"""

def embed_sequence(sequence: str, model_name='esm2_t6_8M_UR50D'):
    try:
        import torch
        import esm
    except ImportError as e:
        raise RuntimeError('Optional ESM-2 support requires: pip install fair-esm') from e
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model.eval()
    batch_converter=alphabet.get_batch_converter(); _,_,tokens=batch_converter([('seq', ''.join(sequence.split()).upper())])
    with torch.no_grad(): out=model(tokens,repr_layers=[model.num_layers],return_contacts=False)
    rep=out['representations'][model.num_layers][0,1:len(sequence.replace(' ',''))+1].mean(0).cpu().numpy()
    return rep.tolist()
