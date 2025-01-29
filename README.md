# parse-who-tb-catalogue-v2
Parse the WHO TB AMR catalogue into GARC for use with gnomonicus

## Output
A parsed, production ready version is available in `NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv`

## Install
Optional virtual environment
```
python -m virtualenv env
source env/bin/activate
```
Install requirements
```
pip install -r requirements.txt
```

## Running
To parse from scratch (not recommended due to processing time):
```bash
python parse-whov2.py --to-garc --parse --filter
```

Because each of the itermediary files is included in the repo, any combination of these options will work (although there is an order to which proccesses which downstream files). e.g `--parse --filter` will produce a re-parsed catalogue and re-filtered interpretation, but `--to-garc --filter` will re-parse coordinates and re-filter, but won't re-parse the catalogue (so re-parsed coordinates won't be used).

Options explained:
* `--to-garc` Parses all rows of the coordinates sheet of the Excel spreadsheet, writing them to `garc_formatted_coordinates.pkl`. This operation takes >55h, so is not recommeded (other than for reproducability) as the output file is included in this repo.
* `--parse` Parses the master sheet of the Excel spreadsheet and the `garc_formatted_coordinates.pkl`. This produces a maximal set of specific rows from the catalogue (`first-pass.csv`). This cannot reliably be used for a few reasons, not excluding the lack of default prediction rules. It's also >52,000 rows so would be unecessarily large.
* `--filter` Filters the specific rows from `first-pass.csv` to remove those covered by default predictions rules (while ensuring the default rules exist). Also ensures that only resistance genes are included, adds minor rules for all group 1/2 specific SNPs, and null rules for all group 1 specific SNPs. There are also a specific set of group 3 (`U`) rules which are covered by default prediction rules, but are added back to increase available evidence for these mutations.


## General strategy
Due to issues found with v1, I decided parsing the ref/alt pairs in the coordinates sheet would be a neat way to wrangle data into a usable format. Many rows in the coordinates sheet specify the same `variant` - showing that these ref/alt pairs are the instances of the given `variant` which have been seen. This leads to a problem; in 1332 `variants` there is no common `variant` between all of the ref/alt pairs. Some of these are vaguely expected, such as non-specific indel `variants`, where all of the parsed rows tend to be more specific. However, there are several cases where I would expect a common `variant`, such as 'premature stop', 'synonymous AA SNP', or 'start lost'. 

To deal with this, we try to catch the rows which a common `variant` is logically unlikely to occur in: 'loss of function', 'amino acid framshift' and 'any non-synonymous AA SNP'. These are then unpacked in a logical way:
* `loss of function` --> [`gene deletion`, `frameshift`, `start SNP`, `premature stop`]
* `amino acid framshift` --> [`framshift at any base of the codon`]
* `any non-synonymous AA SNP` --> [`<gene>@<ref AA><position>*`]

In other cases which logically should have some common `variant`, all variants parsed are added for now. This is not intended to be a long term solution!!

Iterating the master file, we also find `variants` which do not appear within the coordinates sheet. These seem to be restricted to `LoF` and `deletion` mutations though - possibly as a form of 'expert rule'. These are parsed and added separately. 

## Issues (and how this handles them)

### HGVS issues
Some issues arrise from the use of the protein level descriptions of mutations in HGVS. These correspond to pooling of nucleotide level mutations around a common feature to improve the statistical testing for their importance. 

#### AA indel
When inspecting the coordinates sheet for ref/alt pairs, the actual nucleotide indels aren't always as simple as 'insert bases tca at position 3 for a 1_ins_S'. Instead, these insertions can be upstream of the inserted AA, mid codon, but if the post-indel nucleotide sequence is translated, the AA indel behaves as expected.

##### Solution
Unpack the listed AA indels to the instances of nucleotide indels the coordinates sheet details.

#### AA frameshift
For whatever reason, the catalogue includes amino acid level frameshifts. e.g `gid_p.Pro93fs`. This is handled as described above.

##### Solution
Very simply unpack into any frameshift in any base of the codon. e.g `gid@279_fs`, `gid@280_fs`, `gid@281_fs`. Most of these are filtered out due to loss of function rows anyway

### Premature stop codon
I would assume a variant which has a 'premature stop codon' maps to a SNP of `<gene>@*!`. However, the coordinates sheet includes significantly more detailed approach - including a large number of indels which (once AA is translated post-indel), does indeed produce a premature stop codon.
These are flagged as `variants` which do not have a common `variant` between the ref/alt pairs.

##### Solution
Include all mutations parsed from ref/alt pairs. This is not intended to be a long term solution!!

### Start codon lost
For whatever reason (probably convention), all instances of start codons in the catalogue are described as `Met` (`M`). However, TB uses a variety of alternate start codons, so in reality, the start codon is often not `M`.

I would assume this variant would map to a SNP of `<gene>@M1?`, but similar to premature stop, this also includes a variety of upstream indels which lead to a lost start too. 

##### Solution
Include all mutations parsed from ref/alt pairs. This is not intended to be a long term solution!!

### Synonymous AA SNPs
The catalogue returns nucleotide level descriptions of synonymous AA SNPs - being specific about which base was changed. The issue is that there are several cases in which there is no common `variant` between the ref/alt pairs. These often include downstream nucleotide SNPs which can cause the SNP to become non-synonymous at the AA level. This is a known issue raised with the WHO, with a likihood of these rows being removed in future.

##### Solution
Exclude such rows from the catalogue. They are a known bug.

### Catalogue defines intra-gene regions
The catalouge has several (618) ref/alt pairs which do not lie within the gene regions defined by `gumpy`. Most of these are due to very long promoter regions used by the catalogue, which are never picked up by `gumpy` due to an internal limit of 100 bases for a promoter.

##### Solution
Do nothing (for now). Including such rows would require a significant rethink of how `gumpy` assigns promoter regions.



