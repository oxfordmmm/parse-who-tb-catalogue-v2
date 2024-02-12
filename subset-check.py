"""Check specific mutations within WHO v1 against WHO v2
"""
import pandas as pd
import piezo
from collections import defaultdict


def subset(cat: pd.DataFrame) -> dict:
    """Subset a given catalogue to remove default/non-specific rules

    Args:
        cat (pd.DataFrame): Catalogue in

    Returns:
        dict: Mapping of specific rule -> (drug, prediction)
    """
    v1 = defaultdict(dict)
    for i, row in cat.iterrows():
        if "*" in row["MUTATION"] or "?" in row["MUTATION"] or row["MUTATION"].endswith("del"):
            continue
        v1[row["MUTATION"]][row["DRUG"]] = row["PREDICTION"]
    
    return v1

if __name__ == "__main__":
    v1 = subset(pd.read_csv("../who_catalogue_conversion/NC_000962.3_WHO-UCN-GTB-PCI-2021.7_v1.2_GARC1_RFUS.csv"))
    v2 = subset(pd.read_csv("NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv"))
    v2_cat = piezo.ResistanceCatalogue("NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv")
    
    for mutation in sorted(list(v1.keys())):
        pred = v2_cat.predict(mutation)
        same = True
        if pred == "S":
            same = False
            print("No hits in v2 for", mutation, v1[mutation])
        elif len(pred.keys()) != len(v1[mutation].keys()):
            # Different number of drugs
            same = False
            for drug in sorted(list(set(pred.keys()) - set(v1[mutation].keys()))):
                print("Not in v1", mutation, drug, pred[drug])
            for drug in sorted(list(set(v1[mutation].keys()) - set(pred.keys()))):
                print("Not in v2", mutation, drug, v1[mutation][drug])
        else:
            for drug in pred.keys():
                if drug in v1[mutation].keys():
                    if pred[drug] != v1[mutation][drug]:
                        print("Not same prediction for ", mutation, drug, "v1 predicted", v1[mutation][drug], "v2 predicted", pred[drug])
                        same = False
        # if not same:
        #     print()
        


