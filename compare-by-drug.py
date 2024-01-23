"""Find what mutations are different in 2 catalogues for a given drug
"""
import argparse
import pandas as pd

def rename_ordering(ordering: list[str]) -> list[str]:
    """Rename the ordering from filepath to nice names

    Args:
        ordering (list[str]): List of filepaths

    Returns:
        list[str]: List of nice names
    """
    # Assumes filenames are sensible
    new_ordering = []
    for f in ordering:
        if "NC_000962.3_WHO-UCN-GTB-PCI-2021.7_v1.1_GARC1_RFUS.csv" in f:
            new_ordering.append("WHOv1")
        elif "first-pass-filtered.csv" in f:
            new_ordering.append("WHOv2")
        else:
            new_ordering.append(f)
    return new_ordering

def compare(catalogues: dict[str, pd.DataFrame], drug: str) -> None:
    """Find which rules differ for a given drug

    Args:
        catalogues (dict[str, pd.DataFrame]): Mapping of path -> parsed catalogue
        drug (str): 3 letter drug code
    """
    print(drug)
    ordering = sorted(list(catalogues.keys()))
    drug_catalogues = [
        catalogues[catalgoue][catalogues[catalgoue]['DRUG'] == drug][['MUTATION', 'PREDICTION']]
        for catalgoue in ordering
    ]
    cat1 = {
        row['MUTATION']: row['PREDICTION']
        for _, row in drug_catalogues[0].iterrows() 
    }
    cat2 = {
        row['MUTATION']: row['PREDICTION']
        for _, row in drug_catalogues[1].iterrows() 
    }
    both = set()
    just_cat1 = set()
    just_cat2 = set()
    cat1_different_pred = set()
    cat2_different_pred = set()
    for mutation in cat1:
        if cat1[mutation] == "F":
            # v2 doesn't have null rules yet
            continue
        if mutation in cat2.keys():
            # Same mutation in both, check if predictions match
            if cat1[mutation] == cat2[mutation]:
                both.add(mutation)
            else:
                cat1_different_pred.add((mutation, cat1[mutation], cat2[mutation]))
        else:
            just_cat1.add((mutation, cat1[mutation]))

    for mutation in cat2:
        if cat2[mutation] == "F":
            # v2 doesn't have null rules yet
            continue
        if mutation in cat1.keys():
            # Same mutation in both, check if predictions match
            if cat1[mutation] == cat2[mutation]:
                both.add(mutation)
            else:
                cat2_different_pred.add((mutation, cat1[mutation], cat2[mutation]))
        else:
            just_cat2.add((mutation, cat2[mutation]))
    ordering = rename_ordering(ordering)
    print(f"Both catalogues agree on {len(both)} rules!")
    print(f"{ordering[0]} has {len(just_cat1)} more rules")
    print(f"{ordering[1]} has {len(just_cat2)} more rules")
    print(f"{ordering[0]} has different predictions for the same mutation for {len(cat1_different_pred)} rules")
    print(f"{ordering[1]} has different predictions for the same mutation for {len(cat2_different_pred)} rules")
    # print(f"{ordering[0]} extra rules:")
    # for rule in sorted(list(just_cat1)):
    #     print(rule)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogues", action="append", nargs=2, required=True, help="Paths to catalogues to compare")
    parser.add_argument("--drugs", action="append", nargs="+", required=True, help="3 letter drug codes to find results for")
    options = parser.parse_args()

    #Parse the catalogues
    catalogues = {
        catalogue: pd.read_csv(catalogue)
        for catalogue in options.catalogues[0]
    }

    for drug in options.drugs[0]:
        compare(catalogues, drug)
        print()




