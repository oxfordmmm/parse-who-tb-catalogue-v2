"""Parse the WHO catalogue into GARC for use within `gnomonicus`
"""
import argparse
import pickle
import pandas as pd
import numpy as np
import gumpy
from tqdm import tqdm
import asyncio
import tempfile
from collections import Counter, defaultdict
import re
import sys


def build_reference(
    path: str, build: bool = False
) -> tuple[gumpy.Genome, dict[str, gumpy.Gene]]:
    """Build the reference genome object and all genes within it

    Args:
        path (str): Path to genbank file
        build (bool, optional): Whether to build the genome objects from scratch

    Returns:
        gumpy.Genome: Genome object
    """
    if build:
        ref = gumpy.Genome(path)
        ref_genes = {}
        for gene in ref.genes.keys():
            ref_genes[gene] = ref.build_gene(gene)
        pickle.dump(ref_genes, open("reference_genes.pkl", "wb"))
        pickle.dump(ref, open("reference.pkl", "wb"))
    else:
        ref = pickle.load(open("reference.pkl", "rb"))
        ref_genes = pickle.load(open("reference_genes.pkl", "rb"))

    return ref, ref_genes


def background(f):
    def wrapped(*args, **kwargs):
        return asyncio.get_event_loop().run_in_executor(None, f, *args, **kwargs)

    return wrapped


def build_vcf(row: pd.Series) -> gumpy.VCFFile:
    """Parse the ref/alt from the row and build a VCF object for it

    Args:
        row (pd.Series): Row to build a variant from

    Returns:
        gumpy.VCFFile: VCF file object resulting
    """
    pos = row["position"]
    ref = row["reference_nucleotide"]
    alt = row["alternative_nucleotide"]
    vcf = f"""##fileformat=VCFv4.2
##source=minos, version 0.12.5
##fileDate=2023-10-28
##FORMAT=<ID=ALLELE_DP,Number=R,Type=Float,Description="Mean read depth of ref and each allele">
##FORMAT=<ID=COV,Number=R,Type=Integer,Description="Number of reads on ref and alt alleles">
##FORMAT=<ID=COV_TOTAL,Number=1,Type=Integer,Description="Total reads mapped at this site, from gramtools">
##FORMAT=<ID=DP,Number=1,Type=Float,Description="Mean read depth of called allele (ie the ALLELE_DP of the called allele)">
##FORMAT=<ID=FRS,Number=1,Type=Float,Description="Fraction of reads that support the genotype call">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence. Difference in log likelihood of most likely and next most likely genotype">
##minosMeanReadDepth=42.379
##minosReadDepthVariance=937.254
##contig=<ID=NC_000962.3,length=4411532>
##FORMAT=<ID=GT_CONF_PERCENTILE,Number=1,Type=Float,Description="Percentile of GT_CONF">
##FILTER=<ID=MIN_FRS,Description="Minimum FRS of 0.9">
##FILTER=<ID=MIN_DP,Description="Minimum DP of 2">
##FILTER=<ID=MIN_GCP,Description="Minimum GT_CONF_PERCENTILE of 0.5">
##FILTER=<ID=MAX_DP,Description="Maximum DP of 134.22281307415324 (= 3.0 standard deviations from the mean read depth 42.379)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
NC_000962.3	{pos}	.	{ref}	{alt}	.	PASS	.	GT:DP:ALLELE_DP:FRS:COV_TOTAL:COV:GT_CONF:GT_CONF_PERCENTILE	1/1:94:0,94:1.0:94:0,94:590.62:92.99"""
    with tempfile.NamedTemporaryFile("w", delete=False) as f:
        f.write(vcf)
    v = gumpy.VCFFile(str(f.name))
    return v


def parse_ref_alt(
    reference: gumpy.Genome, ref_genes: dict[str, gumpy.Gene], row: pd.Series
) -> list[str]:
    """Use the ref/alt/pos to pull out the variants via VCF

    Args:
        reference (gumpy.Genome): Reference genome
        ref_genes (dict[str, gumpy.Gene]): Reference genes
        row (pd.Series): Row of the catalogue

    Returns:
        list[str]: List of mutations in GARC originating from this row
    """
    ref = row["reference_nucleotide"].lower()
    alt = row["alternative_nucleotide"].lower()
    pos = int(row["position"])

    # Find out which genes this affects
    genes = set()
    for idx, r in enumerate(ref):
        mask = reference.stacked_nucleotide_index == pos + idx
        g = reference.stacked_gene_name[mask]
        for g_ in g:
            genes.add(g_)

    vcf = build_vcf(row)
    sample = reference + vcf
    garc = []
    for gene in sorted(list(genes)):
        if gene:
            ref_gene = ref_genes[gene]
            alt_gene = sample.build_gene(gene)
            diff = ref_gene - alt_gene
            muts = diff.mutations.tolist()
            for mut in muts:
                garc.append(gene + "@" + mut)
    return garc


def shorten_aa(aa: str) -> str:
    """Shorten a 3 letter AA to a single letter AA

    Args:
        aa (str): Long AA

    Returns:
        str: Short AA
    """
    mapping = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
    }
    return mapping[aa.upper()]


def parse_to_garc(
    reference: gumpy.Genome,
    variant: str,
    ref_genes: dict[str, gumpy.Gene],
    row: pd.Series,
) -> list[str]:
    """Parse the given variant to GARC

    Args:
        reference (gumpy.Genome): Reference genome
        variant (str): Variant to parse
        ref_genes (dict[str, gumpy.Gene]): Reference gene objects

    Returns:
        list[str]: Variants in GARC
    """
    garc = []
    # gyrA_LoF
    # This maps to gene deletion, fs, premature stop or start lost
    lof = re.compile(
        r"""
        ([a-zA-Z_0-9.()]+) # Leading gene name
        _LoF # Loss of function
        """,
        re.VERBOSE,
    )
    lof_match = lof.fullmatch(variant)
    if lof_match is not None:
        gene = lof_match.groups()[0]
        # Loss of function maps to a few mutations
        garc.append(gene + "@del_1.0")  # Feature ablation
        garc.append(gene + "@*_fs")  # Frameshift
        if reference.genes[gene]["codes_protein"]:
            garc.append(gene + "@M1?")  # Start lost
            garc.append(gene + "@*!")  # Premature stop
        # print(garc)

    # ***********************************************************************************

    # rrl_n.705A>G | dnaA_c.-46C>G
    n_snp = re.compile(
        r"""
        ([a-zA-Z_0-9.()]+) # Leading gene name
        _[nc]\. # Nucleotide SNP
        (\-?[0-9]+) # Position
        ([ACGT]) # Ref base
        >
        ([ACGT]) # Alt base
        """,
        re.VERBOSE,
    )
    n_snp_match = n_snp.fullmatch(variant)
    if n_snp_match is not None:
        gene, pos, ref, alt = n_snp_match.groups()
        if ref_genes[gene].codes_protein and int(pos) > 0:
            # Catch coding SNPs which for some reason have been allocated nucleotide SNPs here
            garc = parse_ref_alt(reference, ref_genes, row)
        else:
            garc.append(gene + "@" + ref.lower() + pos + alt.lower())

    # ***********************************************************************************

    # Rv1129c_p.Asn5fs
    # I'm assuming this implies a frameshift at any base in the codon
    aa_fs = re.compile(
        r"""
        ([a-zA-Z_0-9.()]+) # Leading gene name
        _p\. # Protein coding
        ([A-Z][a-z][a-z]) # AA 3 letter name
        (\-?[0-9]+) # Position
        fs # Frameshift
        """,
        re.VERBOSE,
    )
    aa_fs_match = aa_fs.fullmatch(variant)
    if aa_fs_match is not None:
        gene, aa, pos = aa_fs_match.groups()
        aa = shorten_aa(aa)
        g = ref_genes[gene]
        r = g.nucleotide_number[g.is_cds][g.codon_number == int(pos)]
        for i in r:
            garc.append(gene + "@" + str(i) + "_fs")
        # print(garc)

    # ***********************************************************************************

    # bacA_p.Met1?
    aa_non_synon = re.compile(
        r"""
        ([a-zA-Z_0-9.()]+) # Leading gene name
        _p\. # Protein coding
        ([A-Z][a-z][a-z]) # AA 3 letter name
        ([0-9]+) # Position
        \? # Non synon
        """,
        re.VERBOSE,
    )
    aa_non_synon_match = aa_non_synon.fullmatch(variant)
    if aa_non_synon_match is not None:
        gene, aa, pos = aa_non_synon_match.groups()
        aa = shorten_aa(aa)
        g = ref_genes[gene]
        r = g.amino_acid_sequence[g.amino_acid_number == int(pos)][0]
        garc.append(gene + "@" + r + "1?")

    # ***********************************************************************************

    aa_snp = re.compile(
        r"""
        ([a-zA-Z_0-9.()]+) # Leading gene name
        _p\. # Protein coding
        ([A-Z][a-z][a-z]) # AA 3 letter name
        ([0-9]+) # Position
        ([A-Z][a-z][a-z]) # AA 3 letter name
        """,
        re.VERBOSE,
    )
    aa_snp_match = aa_snp.fullmatch(variant)
    if aa_snp_match is not None:
        gene, ref_aa, pos, alt_aa = aa_snp_match.groups()
        ref_aa = shorten_aa(ref_aa)
        alt_aa = shorten_aa(alt_aa)
        garc.append(gene + "@" + ref_aa + pos + alt_aa)

    # ***********************************************************************************

    early_stop = re.compile(
        r"""
        ([a-zA-Z_0-9.()]+) # Leading gene name
        _p\. # Protein coding
        ([A-Z][a-z][a-z]) # AA 3 letter name
        ([0-9]+) # Position
        \* # Non synon
        """,
        re.VERBOSE,
    )
    early_stop_match = early_stop.fullmatch(variant)
    if early_stop_match is not None:
        gene, aa, pos = early_stop_match.groups()
        aa = shorten_aa(aa)
        garc.append(gene + "@" + aa + pos + "!")

    if len(garc) == 0:
        # If not already parsed, it's an indel
        #   so avoid nomenclature abiguity by parsing ref/alt
        garc = parse_ref_alt(reference, ref_genes, row)

    return garc


def parse_mutations(scratch: bool = False) -> dict[int, list[str]]:
    """Parse the mutations from the genome coordinates sheet

    Args:
        scratch (bool, optional): Whether to parse from scratch. Takes ~7h
    Returns:
        dict[str, list[str]]: Dict mapping row idx -> list[mutations in garc]
    """
    coordinates = pd.read_excel(
        "WHO-UCN-TB-2023.5-eng.xlsx", sheet_name="Genomic_coordinates"
    )
    if scratch:
        reference, ref_genes = build_reference("NC_000962.3.gbk")
        mutations = {}
        rows = [(idx, row) for idx, row in coordinates.iterrows()]
        for idx, row in tqdm(rows):
            m = parse_to_garc(reference, row["variant"], ref_genes, row)
            mutations[idx] = m
        pickle.dump(mutations, open("garc_formatted2.pkl", "wb"))
    else:
        mutations = pickle.load(open("garc_formatted2.pkl", "rb"))


def parse_confidence_grading(grading: str) -> str:
    """Convert confidence grading to RUS

    Args:
        grading (str): Catalogue grading

    Returns:
        str: Corresponding R, U or S
    """
    conversion = {
        "1) Assoc w R": "R",
        "2) Assoc w R - Interim": "R",
        "3) Uncertain significance": "-",
        "4) Not assoc w R - Interim": "S",
        "5) Not assoc w R": "S",
    }
    return conversion[grading]


def parse_master_file(master_file: pd.DataFrame) -> dict[str, tuple[str, str]]:
    """Parse the master file to pull out the drug/predictions for each

    Args:
        master_file (pd.DataFrame): Master file

    Returns:
        dict[str, tuple[str, str]]: Dict mapping variant -> (drug, pred)
    """
    predictions = defaultdict(list)
    for idx, row in master_file.iterrows():
        drug = row["drug"]
        variant = row["variant"]
        pred = parse_confidence_grading(row["FINAL CONFIDENCE GRADING"])
        predictions[variant].append((drug, pred))
    return predictions


def convert_drug(drug: str) -> str:
    """Convert from the long drug names to 3 letter abbreviations

    Args:
        drug (str): Long named drug

    Returns:
        str: 3 letter drug code
    """
    mapping = {
        "Amikacin": "AMI",
        "Bedaquiline": "BDQ",
        "Capreomycin": "CAP",
        "Clofazimine": "CFZ",
        "Delamanid": "DLM",
        "Ethambutol": "EMB",
        "Ethionamide": "ETH",
        "Isoniazid": "INH",
        "Kanamycin": "KAN",
        "Levofloxacin": "LEV",
        "Linezolid": "LZD",
        "Moxifloxacin": "MXF",
        "Pyrazinamide": "PZA",
        "Rifampicin": "RIF",
        "Streptomycin": "STM",
    }
    return mapping[drug]


def show_problems(
    coordinates: pd.DataFrame, master_file: pd.DataFrame, mutations: list
):
    """Print problem rows (and associated info) to stdout

    Args:
        coordinates (pd.DataFrame): Coordinates sheet
        master_file(pd.DataFrame): Master sheet
        mutations: (list): List of mutations in GARC. Corresponds with rows of the coordinates sheet
    """
    rows = [(idx, row) for idx, row in coordinates.iterrows()]

    common = {}
    all_garc = {}
    not_common = {}
    counting = []
    seen = set()
    for idx, row in rows:
        counting.append(row["variant"])
        if common.get(row["variant"]) is None:
            common[row["variant"]] = set(mutations[idx])
            not_common[row["variant"]] = set(mutations[idx]).difference(
                common[row["variant"]]
            )
            all_garc[row["variant"]] = dict(((idx + 2, mutations[idx]),))
        else:
            common[row["variant"]] = common[row["variant"]].intersection(
                set(mutations[idx])
            )
            not_common[row["variant"]] = (
                not_common[row["variant"]]
                .union(set(mutations[idx]))
                .difference(common[row["variant"]])
            )
            all_garc[row["variant"]][idx + 2] = mutations[idx]
        seen.add(row["variant"])

    for var in common.keys():
        negative = re.compile(r".*-[0-9]+.*")
        negative_match = negative.fullmatch(var)
        if negative_match is not None:
            continue
        pred = parse_confidence_grading(
            master_file[master_file["variant"] == var][
                "FINAL CONFIDENCE GRADING"
            ].values[0]
        )
        if len(common[var]) == 0:
            print(var, all_garc[var], pred, sep="\t")


def catalogue_to_garc(
    rows: list,
    master_file: pd.DataFrame,
    mutations: list,
    reference: gumpy.Genome,
    ref_genes: dict[str, gumpy.Gene],
):
    """Convert the catalogue to useable GARC

    Args:
        rows (list): List of (idx, HGVS) for rows in the coordinates sheet
        master_file (pd.DataFrame): Master file parsed from xlsx
        mutations (list): Mutations parsed from ref/alt pairs
        reference (gumpy.Genome): Reference genome
        ref_genes (dict[str, gumpy.Gene]): Dict mapping gene name -> gumpy Gene object
    """
    seen = set()
    common = {}
    all_garc = {}
    resistance_genes = set()
    resistance_genes_drug = set()
    predictions = parse_master_file(master_file)
    parsed_rows = [
        "GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER\n"
    ]
    COMMON_ALL = "NC_000962.3,WHO-UCN-GTB-PCI-2023.5,1.0,GARC1,RUS,"
    for idx, row in tqdm(rows):
        drug = master_file[master_file["variant"] == row["variant"]]["drug"].values[0]
        pred = parse_confidence_grading(
            master_file[master_file["variant"] == row["variant"]][
                "FINAL CONFIDENCE GRADING"
            ].values[0]
        )
        if pred == "R":
            for m in mutations[idx]:
                g = m.split("@")[0]
                resistance_genes.add(g)
                resistance_genes_drug.add((g, drug))

        if (row["variant"], drug) in seen:
            continue

        # Specific rules which we're ignoring the parsed ref/alt for
        # ***********************************************************************************
        lof = re.compile(
            r"""
            ([a-zA-Z_0-9.()]+) # Leading gene name
            _LoF # Loss of function
            """,
            re.VERBOSE,
        )
        lof_match = lof.fullmatch(row["variant"])
        if lof_match is not None:
            seen.add((row["variant"], drug))
            gene = lof_match.groups()[0]
            garc = []
            # Loss of function maps to a few mutations
            garc.append(gene + "@del_1.0")  # Feature ablation
            garc.append(gene + "@*_fs")  # Frameshift
            if reference.genes[gene]["codes_protein"]:
                first = ref_genes[gene].amino_acid_sequence[0]
                garc.append(gene + "@" + first + "1?")  # Start lost
                garc.append(gene + "@*!")  # Premature stop
            common[row["variant"]] = garc
            if pred == "R":
                for m in garc:
                    g = m.split("@")[0]
                    resistance_genes.add(g)
                    resistance_genes_drug.add((g, drug))
            continue

        # ***********************************************************************************

        aa_fs = re.compile(
            r"""
            ([a-zA-Z_0-9.()]+) # Leading gene name
            _p\. # Protein coding
            ([A-Z][a-z][a-z]) # AA 3 letter name
            (\-?[0-9]+) # Position
            fs # Frameshift
            """,
            re.VERBOSE,
        )
        aa_fs_match = aa_fs.fullmatch(row["variant"])
        if aa_fs_match is not None:
            seen.add((row["variant"], drug))
            gene, aa, pos = aa_fs_match.groups()
            garc = []
            aa = shorten_aa(aa)
            g = ref_genes[gene]
            r = g.nucleotide_number[g.is_cds][g.codon_number == int(pos)]
            for i in r:
                garc.append(gene + "@" + str(i) + "_fs")
            common[row["variant"]] = garc
            if pred == "R":
                for m in garc:
                    g = m.split("@")[0]
                    resistance_genes.add(g)
                    resistance_genes_drug.add((g, drug))
            continue

        # ***********************************************************************************

        aa_non_synon = re.compile(
            r"""
            ([a-zA-Z_0-9.()]+) # Leading gene name
            _p\. # Protein coding
            ([A-Z][a-z][a-z]) # AA 3 letter name
            ([0-9]+) # Position
            \? # Non synon
            """,
            re.VERBOSE,
        )
        aa_non_synon_match = aa_non_synon.fullmatch(row["variant"])
        if aa_non_synon_match is not None:
            seen.add((row["variant"], drug))
            gene, aa, pos = aa_non_synon_match.groups()
            garc = []
            aa = shorten_aa(aa)
            g = ref_genes[gene]
            r = g.amino_acid_sequence[g.amino_acid_number == int(pos)][0]
            garc.append(gene + "@" + r + "1?")
            common[row["variant"]] = garc
            if pred == "R":
                for m in garc:
                    g = m.split("@")[0]
                    resistance_genes.add(g)
                    resistance_genes_drug.add((g, drug))
            continue

        # ***********************************************************************************

        early_stop = re.compile(
            r"""
            ([a-zA-Z_0-9.()]+) # Leading gene name
            _p\. # Protein coding
            ([A-Z][a-z][a-z]) # AA 3 letter name
            ([0-9]+) # Position
            \* # Non synon
            """,
            re.VERBOSE,
        )
        early_stop_match = early_stop.fullmatch(row["variant"])
        if early_stop_match is not None:
            seen.add((row["variant"], drug))
            gene, aa, pos = early_stop_match.groups()
            garc = []
            aa = shorten_aa(aa)
            garc.append(gene + "@" + aa + pos + "!")
            common[row["variant"]] = garc
            if pred == "R":
                for m in garc:
                    g = m.split("@")[0]
                    resistance_genes.add(g)
                    resistance_genes_drug.add((g, drug))
            continue

        # Otherwise, look for a common mutation between all of the parsed rows
        if common.get(row["variant"]) is None:
            common[row["variant"]] = set(mutations[idx])
            all_garc[row["variant"]] = mutations[idx]
        else:
            common[row["variant"]] = common[row["variant"]].intersection(
                set(mutations[idx])
            )
            all_garc[row["variant"]] += mutations[idx]

        for variant in common.keys():
            for drug, pred in predictions[variant]:
                for mutation in common[variant]:
                    g = mutation.split("@")[0]
                    if g in resistance_genes:
                        parsed_rows.append(
                            COMMON_ALL
                            + convert_drug(drug)
                            + ","
                            + mutation
                            + ","
                            + pred
                            + ",{},{},{}\n"
                        )

    with open("default.csv") as f:
        for row in f:
            parsed_rows.append(row.strip())

    with open("first-pass.csv", "w") as f:
        for row in parsed_rows:
            f.write(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--problems",
        action="store_true",
        default=False,
        required=False,
        help="Print information about problem rows to stdout",
    )
    parser.add_argument(
        "--parse",
        action="store_true",
        default=False,
        required=False,
        help="Write the preliminary catalogue, 'first-pass.csv'",
    )
    options = parser.parse_args()
    if not options.problems and not options.parse:
        print("No args given!")
        return

    master_file = pd.read_excel(
        "WHO-UCN-TB-2023.5-eng.xlsx",
        sheet_name="Catalogue_master_file",
        skiprows=[0, 1],
    )
    coordinates = pd.read_excel(
        "WHO-UCN-TB-2023.5-eng.xlsx", sheet_name="Genomic_coordinates"
    )
    try:
        reference, ref_genes = build_reference("NC_000962.3.gbk")
    except:
        reference, ref_genes = build_reference("NC_000962.3.gbk", build=True)
    mutations = pickle.load(open("garc_formatted_all.pkl", "rb"))
    rows = [(idx, row) for idx, row in coordinates.iterrows()]

    if options.problems:
        show_problems(coordinates, master_file, mutations)
    if options.parse:
        catalogue_to_garc(rows, master_file, mutations, reference, ref_genes)


if __name__ == "__main__":
    main()
