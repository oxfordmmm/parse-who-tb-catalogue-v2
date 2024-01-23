"""Compare the performance of the WHOv1 and WHOv2 catalogues
"""
import gnomonicus
import piezo
import pandas as pd
import pickle
import gumpy
from tqdm import tqdm
from collections import defaultdict

def build_reference(path: str) -> tuple[gumpy.Genome, dict[str, gumpy.Gene]]:
    """Build the reference genome object and all genes within it

    Args:
        path (str): Path to genbank file

    Returns:
        gumpy.Genome: Genome object
    """
    # ref = gumpy.Genome(path)
    # ref_genes = {}
    # for gene in ref.genes.keys():
    #     ref_genes[gene] = ref.build_gene(gene)
    # pickle.dump(ref_genes, open("reference_genes.pkl", "wb"))
    # pickle.dump(ref, open("reference.pkl", "wb"))
    ref = pickle.load(open("reference.pkl", "rb"))
    ref_genes = pickle.load(open("reference_genes.pkl", "rb"))
    return ref, ref_genes

def predict_phenotypes(
    mutations: pd.DataFrame,
    catalogue: piezo.ResistanceCatalogue,
    reference_genes: dict[str, gumpy.Gene],
) -> dict[str, str]:
    """Given some mutations, get the effects from the given catalogue, returning a phenotype mapping of drug->pheotype

    Args:
        mutations (pd.DataFrame): Mutations df
        catalogue (piezo.ResistanceCatalogue): Resistance catalgoue
        reference_genes (dict[str, gumpy.Gene]): Reference gene objects

    Returns:
        dict[str, str]: Dictionary mapping drug -> predicted phenotype
    """
    mutations = mutations.rename({'UNIQUEID': 'uniqueid', 'GENE': 'gene', 'MUTATION': 'mutation', 'POSITION': 'position', 'AMINO_ACID_NUMBER': 'amino_acid_number', 'GENOME_INDEX': 'genome_index', 'NUCLEOTIDE_NUMBER': 'nucleotide_number', 'REF': 'ref', 'ALT': 'alt', 'IS_SNP': 'is_snp', 'IS_INDEL': 'is_indel', 'IN_CDS': 'in_cds', 'IN_PROMOTER': 'in_promoter', 'IS_SYNONYMOUS': 'is_synonymous', 'IS_NONSYNONYMOUS': 'is_nonsynonymous', 'IS_HET': 'is_het', 'IS_NULL': 'is_null', 'IS_FILTER_PASS': 'is_filter_pass', 'ELEMENT_TYPE': 'element_type', 'MUTATION_TYPE': 'mutation_type', 'INDEL_LENGTH': 'indel_length', 'INDEL_1': 'indel_1', 'INDEL_2': 'indel_2', 'SITEID': 'siteid', 'NUMBER_NUCLEOTIDE_CHANGES': 'number_nucleotide_changes'}, axis='columns')
    
    #Remove invalid mutations (???)
    to_drop = []
    for idx, row in mutations.iterrows():
        if not reference_genes[row['gene']].valid_variant(row['gene']+"@"+row['mutation']):
            to_drop.append(idx)
    mutations.drop(index=to_drop, inplace=True)
    _, phenotypes, _ = gnomonicus.populateEffects(
        "",
        catalogue,
        mutations,
        reference_genes,
        "",
        False,
        False,
    )
    return {key.split("_")[-1]: phenotypes[key] for key in phenotypes.keys()}, len(to_drop)


def pheno_df_to_dict(phenotypes: pd.DataFrame) -> dict[str, str]:
    """Convert the values given in the phenotypes df to a dict mapping drug -> phenotype

    Args:
        phenotypes (pd.DataFrame): Dataframe of phenotypes for a given sample

    Returns:
        dict[str, str]: Dict mapping drug -> phenotype
    """
    pheno = {}
    for _, row in phenotypes.iterrows():
        drug = row['DRUG']
        p = row['PHENOTYPE']
        if drug in pheno.keys():
            #print("Duplicate drug pheno for ", row['UNIQUEID'])
            if isinstance(p, str) and p in "RFUS"  and "RFUS".index(pheno[drug]) < "RFUS".index(p):
                pheno[drug] = p
        else:
            if isinstance(p, str) and p in "RFUS":
                pheno[drug] = p
    return pheno

def f1_score(tp: int, tn: int, fp: int, fn: int) -> float:
    """Calculate F1 score

    Args:
        tp (int): True positive count
        tn (int): True negative count
        fp (int): False positive count
        fn (int): False negative count

    Returns:
        float: F1 score
    """
    try:
        return tp / (tp + 0.5 * (fp + fn))
    except ZeroDivisionError:
        return 0

def precision(tp: int, tn: int, fp: int, fn: int) -> float:
    """Calculate precision

    Args:
        tp (int): True positive count
        tn (int): True negative count
        fp (int): False positive count
        fn (int): False negative count

    Returns:
        float: precision
    """
    try:
        return tp / (tp + fp)
    except ZeroDivisionError:
        return 0

def accuracy(tp: int, tn: int, fp: int, fn: int) -> float:
    """Calculate accuracy

    Args:
        tp (int): True positive count
        tn (int): True negative count
        fp (int): False positive count
        fn (int): False negative count

    Returns:
        float: accuracy
    """
    try:
        return (tp + tn) / (tp + tn + fp + fn)
    except ZeroDivisionError:
        return 0

def recall(tp: int, tn: int, fp: int, fn: int) -> float:
    """Calculate recall

    Args:
        tp (int): True positive count
        tn (int): True negative count
        fp (int): False positive count
        fn (int): False negative count

    Returns:
        float: recall
    """
    try:
        return tp / (tp + fn)
    except ZeroDivisionError:
        return 0


def compare(phenotypes: list[dict[str, dict[str, str]]], drugs: set[str]) -> None:
    """Compare the performance of the catalogues

    Args:
        phenotypes (list[dict[str, str]]): Phenotypes
    """
    # Positive here is R
    tp_1 = defaultdict(int)
    tn_1 = defaultdict(int)
    fp_1 = defaultdict(int)
    fn_1 = defaultdict(int)

    tp_2 = defaultdict(int)
    tn_2 = defaultdict(int)
    fp_2 = defaultdict(int)
    fn_2 = defaultdict(int)


    for item in phenotypes:
        v1 = item['v1']
        v2 = item['v2']
        real = item['real']

        for drug in real.keys():
            if drug not in drugs:
                continue
            #drugs.add(drug)
            if real[drug] == "R":
                if v1.get(drug,"") == "R":
                    tp_1[drug] += 1
                else:
                    fp_1[drug] += 1
                if v2.get(drug,"") == "R":
                    tp_2[drug] += 1
                else:
                    fp_2[drug] += 1
            else:
                if v1.get(drug,"") != "R":
                    tn_1[drug] += 1
                else:
                    fn_1[drug] += 1
                if v2.get(drug,"") != "R":
                    tn_2[drug] += 1
                else:
                    fn_2[drug] += 1

    for drug in drugs:
        print("************************************************")
        print(drug)
        print("WHO v1:")
        print(f"Precision: {precision(tp_1[drug], tn_1[drug], fp_1[drug], fn_1[drug])}")
        print(f"Recall: {recall(tp_1[drug], tn_1[drug], fp_1[drug], fn_1[drug])}")
        print(f"Accuracy: {accuracy(tp_1[drug], tn_1[drug], fp_1[drug], fn_1[drug])}")
        print(f"F1 score: {f1_score(tp_1[drug], tn_1[drug], fp_1[drug], fn_1[drug])}")
        print()
        print("WHO v2:")
        print(f"Precision: {precision(tp_2[drug], tn_2[drug], fp_2[drug], fn_2[drug])}")
        print(f"Recall: {recall(tp_2[drug], tn_2[drug], fp_2[drug], fn_2[drug])}")
        print(f"Accuracy: {accuracy(tp_2[drug], tn_2[drug], fp_2[drug], fn_2[drug])}")
        print(f"F1 score: {f1_score(tp_2[drug], tn_2[drug], fp_2[drug], fn_2[drug])}")
        print()
        if f1_score(tp_1[drug], tn_1[drug], fp_1[drug], fn_1[drug]) <= f1_score(tp_2[drug], tn_2[drug], fp_2[drug], fn_2[drug]):
            print("WHO v2 is >= WHO v1")
        else:
            print("Worse performance with WHO v2!")
    print("************************************************")
    print("OVERALL")
    tp_1 = sum(tp_1.values())
    tn_1 = sum(tn_1.values())
    fp_1 = sum(fp_1.values())
    fn_1 = sum(fn_1.values())

    tp_2 = sum(tp_2.values())
    tn_2 = sum(tn_2.values())
    fp_2 = sum(fp_2.values())
    fn_2 = sum(fn_2.values())
    print("WHO v1:")
    print(f"Precision: {precision(tp_1, tn_1, fp_1, fn_1)}")
    print(f"Recall: {recall(tp_1, tn_1, fp_1, fn_1)}")
    print(f"Accuracy: {accuracy(tp_1, tn_1, fp_1, fn_1)}")
    print(f"F1 score: {f1_score(tp_1, tn_1, fp_1, fn_1)}")
    print()
    print("WHO v2:")
    print(f"Precision: {precision(tp_2, tn_2, fp_2, fn_2)}")
    print(f"Recall: {recall(tp_2, tn_2, fp_2, fn_2)}")
    print(f"Accuracy: {accuracy(tp_2, tn_2, fp_2, fn_2)}")
    print(f"F1 score: {f1_score(tp_2, tn_2, fp_2, fn_2)}") 



if __name__ == "__main__":
    v1 = piezo.ResistanceCatalogue(
        "NC_000962.3_WHO-UCN-GTB-PCI-2021.7_v1.1_GARC1_RFUS.csv"
    )
    v2 = piezo.ResistanceCatalogue("first-pass-filtered.csv")
    reference, ref_genes = build_reference("NC_000962.3.gbk")

    mutations = pd.read_pickle("MUTATIONS.pkl.gz").reset_index(level=[1,2])
    m_ids = set(mutations.index)
    phenotypes = pd.read_pickle("DST_MEASUREMENTS.pkl.gz").reset_index(level=[1])
    ids = set(phenotypes.index)
    print(f"{len(m_ids)} mutation samples")
    print(f"{len(ids)} phenotype samples")
    print(f"{len(m_ids.intersection(ids))} in both")

    all_phenotypes = []
    skipped = 0
    drugs = sorted(list(set(v1.catalogue.drugs).union(v2.catalogue.drugs)))
    for id_ in tqdm(sorted(list(m_ids.intersection(ids)))):
        muts = mutations[mutations.index == id_]
        predicted_phenotypes_v1, s = predict_phenotypes(muts, v1, ref_genes)
        predicted_phenotypes_v2, s = predict_phenotypes(muts, v2, ref_genes)
        skipped += s
        actual_phenotypes = pheno_df_to_dict(phenotypes[phenotypes.index == id_])
        all_phenotypes.append({
            "v1": predicted_phenotypes_v1,
            "v2": predicted_phenotypes_v2,
            "real": actual_phenotypes
        })
    pickle.dump(all_phenotypes, open("all_phenotypes.pkl", "wb"))
    # all_phenotypes = pickle.load(open("all_phenotypes.pkl", "rb"))
    compare(all_phenotypes, set(drugs))
    print(f"Skippped {skipped} rows due to invalid mutations")
