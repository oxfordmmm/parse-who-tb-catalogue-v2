"""Parse the WHO catalogue into GARC for use within `gnomonicus`
"""
import argparse
import pickle
import pandas as pd
import numpy as np
import gumpy
from tqdm import tqdm
import tempfile
from collections import defaultdict
import re
import json
import piezo


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


def parse_mutations(scratch: bool = False) -> dict[int, list[str]]:
    """Parse the mutations from the genome coordinates sheet

    Args:
        scratch (bool, optional): Whether to parse from scratch. Takes ~7h
    Returns:
        dict[str, list[str]]: Dict mapping row idx -> list[mutations in garc]
    """
    coordinates = pd.read_excel(
        "WHO-UCN-TB-2023.7-eng.xlsx", sheet_name="Genomic_coordinates"
    )
    if scratch:
        reference, ref_genes = build_reference("NC_000962.3.gbk")
        mutations = {}
        rows = [(idx, row) for idx, row in coordinates.iterrows()]
        for idx, row in tqdm(rows):
            m = parse_ref_alt(reference, ref_genes, row)
            mutations[idx] = m
        pickle.dump(mutations, open("garc_formatted_coordinates.pkl", "wb"))
    else:
        mutations = pickle.load(open("garc_formatted_coordinates.pkl", "rb"))


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
        "3) Uncertain significance": "U",
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


def check(a, b):
    if np.isnan(a) or np.isnan(b):
        return
    # if a != b:
    #     print(a, b)
    assert a >= b, f"{a} != {b}"


def safe_solo(solo) -> int | None:
    """Take a value from the master sheet and conver to None if applicable

    Args:
        solo (Any): Solo from master sheet

    Returns:
        int | None: Integer if int, None if nan
    """
    if isinstance(solo, str):
        return solo
    if np.isnan(solo):
        return None
    return int(solo)


def get_evidences(master_file: pd.DataFrame) -> tuple[dict, dict]:
    """Parse the evidence fields out of the master sheet.

    Args:
        master_file (pd.DataFrame): Master file sheet

    Returns:
        tuple[dict, dict]: Tuple of (dict mapping (HGVS variant, drug) -> evidence JSON, dict mapping (HGVS variant, drug) -> other JSON)
    """
    evidences = {}
    others = {}
    for _, row in master_file.iterrows():
        variant = row["variant"]
        drug = convert_drug(row["drug"])
        if evidences.get((variant, drug)) is not None:
            pred = parse_confidence_grading(row["FINAL CONFIDENCE GRADING"])
            print(
                f">1 row for {variant} -> {drug}. This one is {pred}. Last is {parse_confidence_grading(others.get((variant, drug))['FINAL_CONFIDENCE_GRADING'])}"
            )

        evidences[(variant, drug)] = {
            "Observed_samples": {
                "Present_SOLO_SR": safe_solo(row["Present_SOLO_SR"]),
                "Present_SOLO_R": safe_solo(row["Present_SOLO_R"]),
                "Present_SOLO_S": safe_solo(row["Present_SOLO_S"]),
            },
            "Additional grading criteria applied": safe_solo(
                row["Additional grading criteria applied"]
            ),
            "FINAL CONFIDENCE GRADING": safe_solo(row["FINAL CONFIDENCE GRADING"]),
            "INITIAL CONFIDENCE GRADING": safe_solo(row["INITIAL CONFIDENCE GRADING"]),
            "WHO HGVS": row["variant"],
        }
    return evidences, others


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
    evidences, others = get_evidences(master_file)
    seen = set()
    common = {}
    all_garc = {}
    resistance_genes = set()
    resistance_genes_drug = set()
    predictions = parse_master_file(master_file)
    parsed_rows = []
    COMMON_ALL = ["NC_000962.3", "WHO-UCN-GTB-PCI-2023.5", "1.0", "GARC1", "RUS"]
    skipped = []
    for idx, row in tqdm(rows):
        if len(master_file[master_file["variant"] == row["variant"]]) == 0:
            # Variant isn't in the master file ??????
            skipped.append(row["variant"])
            continue
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
            if len(common[variant]) == 0:
                # If we've hit a problem row, exclude non-synon SNPs (due to known issue),
                # then just take coordinate rows for now
                nc_snp = re.compile(
                    r"""
                    ([a-zA-Z_0-9.()]+) # Leading gene name
                    _c\.
                    ([0-9]+) # Position
                    ([ACGT]) # Ref base
                    >
                    ([ACGT]) # Alt base
                """,
                    re.VERBOSE,
                )
                nc_snp_match = nc_snp.fullmatch(variant)
                if nc_snp_match is not None:
                    continue
                common[variant] = sorted(list(set(all_garc[variant])))

            # Use multi-mutations where applicable
            # (rules we manually unpack are lists, ref/alt rules are sets)
            if isinstance(common[variant], set):
                # Check for large dels (especially ones parsed from ref/alt)
                # These will occur with specific dels and shouldn't be joined
                large_del = re.compile(
                    r"""
                    ([a-zA-Z_0-9.()]+) # Leading gene name
                    @del_0\.([0-9][0-9]) # Any large del
                    """,
                    re.VERBOSE,
                )
                safe_mutations = []
                for mutation in common[variant]:
                    if "Rv2042c" in mutation:
                        # These are introduced due to pncA_p.Ter187Trpext*?
                        # so we shouldn't consider Rv2042c as a resistance gene
                        continue
                    large_del_match = large_del.fullmatch(mutation)
                    if large_del_match is not None:
                        # We have a large del so we only really care about this
                        safe_mutations = [mutation]
                        break
                    safe_mutations.append(mutation)

                mutation = "&".join(sorted(list(safe_mutations)))
                parsed_rows.append(
                    [
                        *COMMON_ALL,
                        convert_drug(drug),
                        mutation,
                        pred,
                        json.dumps({}),
                        json.dumps(evidences.get((variant, convert_drug(drug)), {})),
                        json.dumps(others.get((variant, convert_drug(drug)), {})),
                    ]
                )
            else:
                lof = re.compile(
                    r"""
                    ([a-zA-Z_0-9.()]+) # Leading gene name
                    _LoF # Loss of function
                    """,
                    re.VERBOSE,
                )
                lof_match = lof.fullmatch(variant)
                for mutation in common[variant]:
                    parsed_rows.append(
                        [
                            *COMMON_ALL,
                            convert_drug(drug),
                            mutation,
                            pred,
                            json.dumps({}),
                            json.dumps(
                                evidences.get((variant, convert_drug(drug)), {})
                            ),
                            json.dumps(others.get((variant, convert_drug(drug)), {})),
                        ]
                    )

    expert = pd.read_csv("expert-rules.csv")
    report_expert = pd.read_csv("report-expert-rules.csv")
    cat = pd.DataFrame(
        parsed_rows,
        columns=[
            "GENBANK_REFERENCE",
            "CATALOGUE_NAME",
            "CATALOGUE_VERSION",
            "CATALOGUE_GRAMMAR",
            "PREDICTION_VALUES",
            "DRUG",
            "MUTATION",
            "PREDICTION",
            "SOURCE",
            "EVIDENCE",
            "OTHER",
        ],
    )
    cat = pd.concat([cat, expert, report_expert])
    cat.to_csv("first-pass.csv", index=False)
    return skipped


def fetch_group_3_to_keep() -> set[tuple[str, str]]:
    """Parse the group 3 rules to keep (hand picked by Tim Peto)

    Note that these have 0 impact on prediction as they're all U's,
    but improve evidence when these rules are hit

    Returns:
        set[tuple[str, str]]: Set of (WHO_HGVS, drug) tuples to keep
    """
    group_3 = pd.read_csv("group-3-to-keep.csv")
    to_keep = set()
    for _, row in group_3.iterrows():
        to_keep.add((row["variant"], convert_drug(row["drug"])))
    return to_keep


def filter(reference_genes: dict[str, gumpy.Gene]):
    """Remove superfluous rows which are covered by default or broad rules
    Almost exactly the same proceedure as with v1. Also adds default rules

    Args:
        reference_genes (dict[str, gumpy.Gene]): Dict mapping gene name -> gene object
    """
    catalogue = pd.read_csv("first-pass.csv")
    resistanceGenes = {("mmpL5", "BDQ"), ("mmpL5", "CFZ")}
    seen = set()
    lof_genes = {}

    # Find all of the genes which confer resistance to a given drug
    for i, row in catalogue.iterrows():
        prediction = row["PREDICTION"]
        mutation = row["MUTATION"]
        drug = row["DRUG"]
        if "LoF" in json.loads(row["EVIDENCE"]).get("WHO HGVS", []):
            lof_genes[(mutation.split("@")[0], drug)] = prediction
        if prediction == "R":
            if "&" in mutation:
                if mutation[0] == "^":
                    mutation = mutation[1::]
                for m in mutation.split("&"):
                    resistanceGenes.add((m.split("@")[0], drug))
            else:
                resistanceGenes.add((mutation.split("@")[0], drug))
    to_keep = fetch_group_3_to_keep()
    fixed = {col: [] for col in catalogue}
    for i, row in catalogue.iterrows():
        toDelete = False
        prediction = row["PREDICTION"]
        mutation = row["MUTATION"]
        if (mutation, row["DRUG"], prediction) in seen:
            # Avoid duplicate rows
            print(
                f"Removing {row['MUTATION']}:{row['DRUG']}:{row['PREDICTION']} as it already exists!"
            )
            toDelete = True
        if "&" in mutation:
            if mutation[0] == "^":
                toDelete = False
            else:
                toDelete = False
                for mut in mutation.split("&"):
                    if (mut.split("@")[0], row["DRUG"]) not in resistanceGenes:
                        print(
                            f"Removing {mutation} --> {row['DRUG']} : {(mut.split('@')[0], row['DRUG'])} as it's not a resistance gene for {mut}"
                        )
                        toDelete = True
        if (
            "^" not in mutation
            and (mutation.split("@")[0], row["DRUG"]) not in resistanceGenes
        ):
            toDelete = True
            print(
                f"Removing {row['MUTATION']}:{row['DRUG']}:{row['PREDICTION']} as it is not a resistance gene"
            )

        elif prediction == "U":
            if (
                (json.loads(row["EVIDENCE"]).get("WHO HGVS"), row["DRUG"])
                not in to_keep
                or "dup" in json.loads(row["EVIDENCE"]).get("WHO HGVS")
                or "ext" in json.loads(row["EVIDENCE"]).get("WHO HGVS")
            ):
                indel = re.compile(
                    r"""
                                ([a-zA-Z_0-9]+@) #Leading gene name
                                (
                                    (-?[0-9]+_((ins)|(del))_[acgotxz]*) #indel
                                )
                                """,
                    re.VERBOSE,
                )
                if indel.fullmatch(mutation):
                    # Matched an indel generic so delete
                    toDelete = True
                    print(
                        f"Removing {row['MUTATION']}:{row['DRUG']}:{row['PREDICTION']} as it matches *_indel-->U"
                    )
                nonsynon = re.compile(
                    r"""
                                    ([a-zA-Z_0-9]+@) #Leading gene name
                                    (([!ACDEFGHIKLMNOPQRSTVWXYZacgotxz])-?[0-9]+([!ACDEFGHIKLMNOPQRSTVWXYZacgotxz])) #SNP
                                    """,
                    re.VERBOSE,
                )
                if nonsynon.fullmatch(mutation):
                    _, _, base1, base2 = nonsynon.fullmatch(mutation).groups()
                    if base1 != base2:
                        # This matches the gene@*? or gene@-*? so delete
                        toDelete = True
                        print(
                            f"Removing {row['MUTATION']}:{row['DRUG']}:{row['PREDICTION']} as it matches *?-->U or -*?-->U"
                        )
            if mutation.split("@")[0] == "mmpL5":
                toDelete = True

        elif prediction == "S":
            # Checking for gene@*=
            # This has now become gene@A*A(&gene@<nucleotide><pos><nucleotide>){1,3}
            synon = re.compile(
                r"""
                                ([a-zA-Z_0-9]+@) #Leading gene name
                                (([!ACDEFGHIKLMNOPQRSTVWXYZ])[0-9]+([!ACDEFGHIKLMNOPQRSTVWXYZ])) #SNP
                                (& #And the nucleotide(s) causing this
                                ([a-zA-Z_0-9]+@) #Gene
                                ([a-z][0-9]+[a-z]))+ #Nucleotide SNP
                                """,
                re.VERBOSE,
            )
            if synon.fullmatch(mutation):
                matches = synon.fullmatch(mutation).groups()
                base1 = matches[2]
                base2 = matches[3]
                if base1 == base2:
                    # Matches the synonymous mutation so delete
                    toDelete = True
                    print(
                        f"Removing {row['MUTATION']}:{row['DRUG']}:{row['PREDICTION']} as it matches *=-->S"
                    )

        if lof_genes.get((mutation.split("@")[0], row["DRUG"])) is not None:
            # This matches a default rule for LoF
            # Fist check this predicts <= LoF for this gene
            if "RUS".index(
                lof_genes[(mutation.split("@")[0], row["DRUG"])]
            ) >= "RUS".index(prediction):
                fs = re.compile(
                    r"""
                    ([a-zA-Z_0-9]+@) #Leading gene name
                    [0-9]+ #pos
                    _fs
                    """,
                    re.VERBOSE,
                )
                fs_match = fs.fullmatch(mutation)
                if fs_match is not None:
                    # This matches the default rule so remove it
                    toDelete = True
                    print(
                        f"Removing {mutation}:{row['DRUG']}:{prediction} as it matches {mutation.split('@')[0]}_LoF:{lof_genes.get((mutation.split('@')[0], row['DRUG']))}"
                    )

        if not toDelete:
            # We want to keep this one, so add to fixed
            for col in row.axes[0]:
                fixed[col].append(row[col])

        seen.add((mutation, row["DRUG"], prediction))

    # Add default rules for all resistance genes
    for gene, drug in sorted(list(resistanceGenes)):
        if gene == "mmpL5":
            pred = "S"
        else:
            pred = "U"
        defaults = [
            (gene + "@*?", pred),
            (gene + "@-*?", pred),
            (gene + "@*_indel", pred),
            (gene + "@-*_indel", pred),
            (gene + "@del_0.0", pred),
        ]
        if reference_genes[gene].codes_protein:
            defaults.append((gene + "@*=", "S"))
        for mut, pred in defaults:
            fixed["GENBANK_REFERENCE"].append("NC_000962.3")
            fixed["CATALOGUE_NAME"].append("WHO-UCN-GTB-PCI-2023.5")
            fixed["CATALOGUE_VERSION"].append("1.0")
            fixed["CATALOGUE_GRAMMAR"].append("GARC1")
            fixed["PREDICTION_VALUES"].append("RUS")
            fixed["DRUG"].append(drug)
            fixed["MUTATION"].append(mut)
            fixed["PREDICTION"].append(pred)
            fixed["SOURCE"].append("{}")
            fixed["EVIDENCE"].append('{"reporting_rule": "True"}')
            fixed["OTHER"].append("{}")

    catalogue = pd.DataFrame(fixed)

    catalogue.to_csv("first-pass-filtered.csv", index=False)
    add_nulls_and_minors()


def filter_evidence(r_mutations: list[dict]) -> dict:
    if len(r_mutations) == 1:
        return r_mutations[0]
    most_significant_mut = None
    most_significant_total = 0
    for r_mut, evidence in r_mutations:
        ev = evidence.get("Observed_samples", {})
        solo_r = ev.get("Present_SOLO_R", "-")
        solo_sr = ev.get("Present_SOLO_SR", "-")
        if solo_r == "-" or solo_sr == "-":
            # No proper evidence so skip
            continue
        if solo_r is None or solo_sr is None:
            # No proper evidence so skip
            continue
        if (solo_r + solo_sr) >= most_significant_total:
            # most_significant_ev = percent
            most_significant_mut = (r_mut, evidence)
            most_significant_total = solo_r + solo_sr

    if most_significant_mut is None:
        # This only happens with expert rules not otherwise in catalogue
        # So just take the expert rule (there's 1 mutation here)
        most_significant_mut = r_mutations[0]
        # print("@@@", r_mutations)

    return most_significant_mut


def add_nulls_and_minors():
    path = "first-pass-filtered.csv"
    original = piezo.ResistanceCatalogue(path)
    original_df = pd.read_csv(path)
    rules = original.catalogue.rules

    snps = rules[rules["MUTATION_TYPE"] == "SNP"]
    r_snps = snps[snps["PREDICTION"] == "R"]
    specific_r_snps = r_snps[r_snps["POSITION"] != "*"]
    indels = rules[rules["MUTATION_TYPE"] == "INDEL"]
    r_indels = indels[indels["PREDICTION"] == "R"]
    f_rules = defaultdict(list)
    minor_rules = defaultdict(list)

    # Add minor rules for all R SNPs (including wildcards)
    for _, row in r_snps.iterrows():
        if str(row["MUTATION"][-1]).islower():
            # Nucleotide mutation
            null_mutation = row["MUTATION"][:-1] + "x"
        else:
            null_mutation = row["MUTATION"][:-1] + "X"

        if row["EVIDENCE"].get("FINAL CONFIDENCE GRADING", "") == "1) Assoc w R":
            
            minor_rules[(row["GENE"], row["MUTATION"] + ":3", row["DRUG"])].append(
                (row["GENE"] + "@" + row["MUTATION"], row["EVIDENCE"])
            )
        if (
            row["EVIDENCE"].get("FINAL CONFIDENCE GRADING", "")
            == "2) Assoc w R - Interim"
        ):
            minor_rules[(row["GENE"], row["MUTATION"] + ":3", row["DRUG"])].append(
                (row["GENE"] + "@" + row["MUTATION"], row["EVIDENCE"])
            )

    # Add minor rules for all R indels (including wildcards)
    for _, row in r_indels.iterrows():
        minor_rules[(row["GENE"], row["MUTATION"] + ":3", row["DRUG"])].append(
            (row["GENE"] + "@" + row["MUTATION"], row["EVIDENCE"])
        )

    # Only add null rules for specific R SNPs
    for _, row in specific_r_snps.iterrows():
        if str(row["MUTATION"][-1]).islower():
            # Nucleotide mutation
            null_mutation = row["MUTATION"][:-1] + "x"
        else:
            null_mutation = row["MUTATION"][:-1] + "X"

        if row["EVIDENCE"].get("FINAL CONFIDENCE GRADING", "") == "1) Assoc w R":
            # Only F specifics
            f_rules[(row["GENE"], null_mutation, row["DRUG"])].append(
                (row["GENE"] + "@" + row["MUTATION"], row["EVIDENCE"])
            )

    f_df_dict = defaultdict(list)

    print(f"Adding {len(f_rules.keys())} null rules")
    for g, m, d in f_rules.keys():
        f_df_dict["GENBANK_REFERENCE"].append("NC_000962.3")
        f_df_dict["CATALOGUE_NAME"].append("WHO-UCN-GTB-PCI-2023.5")
        f_df_dict["CATALOGUE_VERSION"].append(2.0)
        f_df_dict["CATALOGUE_GRAMMAR"].append("GARC1")
        f_df_dict["PREDICTION_VALUES"].append("RFUS")
        f_df_dict["DRUG"].append(d)
        f_df_dict["MUTATION"].append(g + "@" + m)
        f_df_dict["PREDICTION"].append("F")
        f_df_dict["SOURCE"].append({})
        significant_mut, significant_ev = filter_evidence(f_rules[(g, m, d)])
        significant_mut = significant_mut.split("@")[1]
        f_df_dict["EVIDENCE"].append(
            json.dumps({"Resistant prediction": significant_mut, **significant_ev})
        )
        f_df_dict["OTHER"].append({})

    print(f"Adding {len(minor_rules.keys())} minor rules")
    for g, m, d in minor_rules.keys():
        f_df_dict["GENBANK_REFERENCE"].append("NC_000962.3")
        f_df_dict["CATALOGUE_NAME"].append("WHO-UCN-GTB-PCI-2023.5")
        f_df_dict["CATALOGUE_VERSION"].append(2.0)
        f_df_dict["CATALOGUE_GRAMMAR"].append("GARC1")
        f_df_dict["PREDICTION_VALUES"].append("RFUS")
        f_df_dict["DRUG"].append(d)
        f_df_dict["MUTATION"].append(g + "@" + m)
        f_df_dict["PREDICTION"].append("R")
        f_df_dict["SOURCE"].append({})
        significant_mut, significant_ev = filter_evidence(minor_rules[(g, m, d)])
        significant_mut = significant_mut.split("@")[1]
        f_df_dict["EVIDENCE"].append(json.dumps(significant_ev))
        f_df_dict["OTHER"].append({})

    f_df = pd.DataFrame.from_dict(f_df_dict)

    new_df = pd.concat([original_df, f_df])
    new_df["CATALOGUE_VERSION"] = 2.0
    new_df["PREDICTION_VALUES"] = "RFUS"
    new_df.reset_index(drop=True, inplace=True)
    new_df.to_csv(
        "NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv", index=False
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--to-garc",
        action="store_true",
        default=False,
        required=False,
        help="Convert the genome coordinates sheet to GARC. Takes ~55h.",
    )
    parser.add_argument(
        "--parse",
        action="store_true",
        default=False,
        required=False,
        help="Write the preliminary catalogue, 'first-pass.csv'. Takes ~10m.",
    )
    parser.add_argument(
        "--filter",
        action="store_true",
        default=False,
        required=False,
        help="Filter out rows covered by general rules.",
    )
    options = parser.parse_args()
    if not options.to_garc and not options.parse and not options.filter:
        print("No args given!")
        return

    if options.to_garc:
        # Double check that the user knows that a) this is destructive and b) takes a long time
        check = input(
            "The to-garc operation will take ~55h and overwrite the existing garc_formatted_coordinates.pkl, continue? [y/n]"
        )
        if check == "y":
            parse_mutations(scratch=True)
        else:
            print("Aborting!")
            return

    # Get the reference genome and genes. Trying from pkl first for speed.
    try:
        reference, ref_genes = build_reference("NC_000962.3.gbk")
    except:  # noqa: E722
        reference, ref_genes = build_reference("NC_000962.3.gbk", build=True)

    if options.parse:
        master_file = pd.read_excel(
            "WHO-UCN-TB-2023.7-eng.xlsx",
            sheet_name="Catalogue_master_file",
            skiprows=[0, 1],
        )
        coordinates = pd.read_excel(
            "WHO-UCN-TB-2023.7-eng.xlsx", sheet_name="Genomic_coordinates"
        )

        mutations = pickle.load(open("garc_formatted_coordinates.pkl", "rb"))
        rows = [(idx, row) for idx, row in coordinates.iterrows()]
        catalogue_to_garc(rows, master_file, mutations, reference, ref_genes)

    if options.filter:
        filter(ref_genes)


if __name__ == "__main__":
    main()
