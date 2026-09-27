from preprocessing.structure_processing import lysine_coordinates

def summarize_structure(path):
    lys=lysine_coordinates(path)
    return {"structure":str(path),"lysine_count":len(lys),"chains":sorted({x["chain"] for x in lys})}
