import pandas as pd

def merge_label_to_fm(fm,label):
    
    fm = fm.copy()
    fm["label"] = label

    fm = fm[fm["label"].isin(["R", "S"])]
    fm["label"] = fm["label"].map({"R": 1, "S": 0})

    return fm