
import pickle
import gumpy
import tempfile
import pandas as pd
from tqdm import tqdm

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


def build_vcf(ref: str, alt: str, pos: int) -> gumpy.VCFFile:
    """Parse the ref/alt from the row and build a VCF object for it

    Args:
        row (pd.Series): Row to build a variant from

    Returns:
        gumpy.VCFFile: VCF file object resulting
    """
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
    reference: gumpy.Genome,
    ref_genes: dict[str, gumpy.Gene],
    ref: str,
    alt: str,
    pos: int,
) -> list[str]:
    """Use the ref/alt/pos to pull out the variants via VCF

    Args:
        reference (gumpy.Genome): Reference genome
        ref_genes (dict[str, gumpy.Gene]): Reference genes
        row (pd.Series): Row of the catalogue

    Returns:
        list[str]: List of mutations in GARC originating from this row
    """

    # Find out which genes this affects
    genes = set()
    for idx, r in enumerate(ref):
        mask = reference.stacked_nucleotide_index == pos + idx
        g = reference.stacked_gene_name[mask]
        for g_ in g:
            genes.add(g_)

    vcf = build_vcf(ref, alt, pos)
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

if __name__ == "__main__":
    reference, ref_genes = build_reference("NC_000962.3.gbk")
    
    old_cat = pd.read_excel("WHO-UCN-TB-2023.5-eng.xlsx", sheet_name="Genomic_coordinates")
    new_cat = pd.read_excel("WHO-UCN-TB-2023.6-eng.xlsx", sheet_name="Genomic_coordinates")

    old = set([tuple(row.tolist()) for _, row in old_cat.iterrows()])
    new = set([tuple(row.tolist()) for _, row in new_cat.iterrows()])

    to_add = new - old
    to_skip = old - new
    print(len(to_add))
    x = [row[0] for row in to_add]
    pickle.dump(x, open("test.pkl", "wb"))

    old_map = {tuple(row.tolist()): idx for idx, row in old_cat.iterrows()}
    new_rows = [(idx, tuple(row.tolist())) for idx, row in new_cat.iterrows()]

    # This is a mapping of old row number --> GARC for that row
    old_garc = pickle.load(open("garc_formatted_all.pkl", "rb"))
    # We want to skip over ones not in the new one
    to_skip_idx = [old_map[row] for row in sorted(to_skip)]
    for idx in to_skip_idx:
        del old_garc[idx]

    updated_rows = pickle.load(open("updated-rows.pkl", "rb"))

    new_garc = {}

    for idx, row in tqdm(new_rows):
        if old_map.get(row) is None:
            # This is a new row so pull from the new rows
            new_garc[idx] = updated_rows[row]
        else:
            # Already existed so skip it
            new_garc[idx] = old_garc[old_map[row]]
    
    # pickle.dump(new_garc, open("garc_formatted_8_feb_2024.pkl", "wb"))

    # parsed = {}
    # for variant, chrom, pos, ref, alt in tqdm(sorted(list(to_add))):
    #     g = parse_ref_alt(reference, ref_genes, ref, alt, int(pos))
    #     parsed[(variant, chrom, pos, ref, alt)] = g
    # pickle.dump(parsed, open("updated-rows.pkl", "wb"))

