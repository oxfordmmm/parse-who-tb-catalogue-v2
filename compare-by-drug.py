"""Find what mutations are different in 2 catalogues for a given drug
"""
import argparse
import pandas as pd

def classify_rule(row: pd.core.series.Series) -> str:
    rule = row.gene_mutation
    if '*' in rule:
        return 'default'
    elif '&' in rule:
        return 'combination'
    elif '_fs' in rule:
        return 'frameshift'
    elif '_ins' in rule:
        return 'insertion'
    elif '_del' in rule:
        return 'deletion'
    elif 'del' in rule:
        return 'deletion'
    elif '-' in rule:
        return 'promoter'
    else:
        return 'snp'

def create_df(both: set[str], just_cat1: set[str], just_cat2: set[str]) -> pd.DataFrame:
    """Create a DataFrame with the rules and their catalogue membership for further analysis
    
    Args:
        both (set[str]): Rules that are in both catalogues
        just_cat1 (set[str]): Rules that are only in catalogue 1
        just_cat2 (set[str]): Rules that are only in catalogue 2

    Returns:
        pandas.DataFrame: Dataframe of rules and their catalogue membership
    """
    df_both = pd.DataFrame(list(both), columns=['gene_mutation', 'phenotype'])
    df_both['membership'] = 'both'
    df_both['cat1'] = True
    df_both['cat2'] = True

    df_just_cat1 = pd.DataFrame(list(just_cat1), columns=['gene_mutation', 'phenotype'])
    df_just_cat1['membership'] = 'cat1'
    df_just_cat1['cat1'] = True
    df_just_cat1['cat2'] = False

    df_just_cat2 = pd.DataFrame(list(just_cat2), columns=['gene_mutation', 'phenotype'])
    df_just_cat2['membership'] = 'cat2'
    df_just_cat2['cat1'] = False
    df_just_cat2['cat2'] = True

    df = pd.concat([df_both, df_just_cat1, df_just_cat2])

    df['rule-type'] = df.apply(classify_rule, axis=1)

    return df

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

    Returns:
        pandas.DataFrame: Dataframe of rules and their catalogue membership
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
                both.add((mutation, cat1[mutation]))
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
                both.add((mutation, cat1[mutation]))
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

    if options.write_csv:
        df = create_df(both, just_cat1, just_cat2)
        df['drug'] = drug
        df = df[['drug', 'membership','cat1','cat2','rule-type', 'phenotype', 'gene_mutation']]
        return df
    else:
        return None
    # print(f"{ordering[0]} extra rules:")
    # for rule in sorted(list(just_cat1)):
    #     print(rule)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogues", action="append", nargs=2, required=True, help="Paths to catalogues to compare")
    parser.add_argument("--drugs", action="append", nargs="+", required=True, help="3 letter drug codes to find results for")
    parser.add_argument("--write-csv", action='store_true', help="write out the rules and their catalogue membership to a CSV file for further analysis")
    options = parser.parse_args()

    #Parse the catalogues
    catalogues = {
        catalogue: pd.read_csv(catalogue)
        for catalogue in options.catalogues[0]
    }

    if options.drugs[0][0] == 'ALL':
        dfs = []
        drugs = set()
        for cat in catalogues:
            for drug in catalogues[cat].DRUG.unique():
                drugs.add(drug)
        drugs = sorted(list(drugs))
    else:
        drugs = options.drugs[0]


    for drug in drugs:
        if options.drugs[0][0] != 'ALL':
            df = compare(catalogues, drug)
            df.to_csv(f"{drug}.csv", index=False)
        else:
            dfs.append(compare(catalogues, drug))

        print()

    if options.drugs[0][0] == 'ALL':
        df = pd.concat(dfs)
        df.to_csv(f"ALL.csv", index=False)




