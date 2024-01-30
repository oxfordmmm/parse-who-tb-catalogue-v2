# parse-who-tb-catalogue-v2
Parse the WHO TB AMR catalogue into GARC for use with gnomonicus

## Initial version
A possibly usable first draft is available in `NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv`

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
This is still a WIP, so less than ideal. 
1. Parse every row of the `coordinates` sheet with gumpy. This takes ~3 days. Code is within `parse.py` for this, but not currently available as a CLI option
2. Parse the master file and build a CSV of rows which don't exist in the coordinates sheet: `python parse.py --t`. This produces `expert-rules.csv`
3. Parse a very verbose and literal translation of the master sheet: `python parse.py --parse`. This produces `first-pass.csv`
4. Filter out unhelpful rows and add default rules: `python parse.py --filter`. This produces `first-pass-filtered.csv` 


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


## Known issues
* No minor population rules yet


