# DAYWWYD (ver. B)
## Do Anything You Want With Your Data (Better than ver. A)
*Expect a new version to be released: L (Literally). Maybe.*
> *This is the repo for the homework of the BI Python 2023 course*

**DAYWWYD is set of scripts for different bioinformatics operations.**

It consit of two main scrits:
- `DAYWWYD.A.py`
- `bio_files_processor.py`

## Installation
`git clone git@github.com:KirPetrikov/DAYWWYD.git`

## Requirements and dependencies
*All is included in the package. Please don't remove `src` folder or any of its content*

`DAYWWYD.A` requires:
- `na_seq_tool.py`
- `prot_seq_tool.py`
- `fastq_filter.py`
- `parse_fastq.py`
- `file_names.py`

`bio_files_processor` requires:
- `parse_gbk.py`
- `file_names.py`

Service modules should be located in the `./src` folder.

## `DAYWWYD.A.py`

It includes three functions:
- `process_na` operations on nucleic acids sequences
- `process_prot` operations on proteins (amino acids sequenses)
- `filter_fastq` filtering sequences from fastq files

### `process_na(operation, seqs)`

Accepts two variables:
- `operation` defines defines the operation to be applied
- `seqs` list of strings with sequences to process
Returns list of perocessed sequences, the letters case is preserved.

**Valid characters** for sequences: `A, T, G, C, U, a, t, g, c, u`. Chimeric DNA-RNA sequences, i.e. including both `T, t`, and `U, u` aren't allowed.

**Operations**:
- `trans` returns the transcribed sequence. If the original sequence contained `U, u` then they will be replaced by `T, t` as for reverse transcription
- `rev` returns the inverted sequence, from backward to forward
- `comp` returns the complementary sequence
- `revcomp` returns the inverted complementary sequence

**Examples**
```python
process_na('trans', ['ATG']) # [AUG]
process_na('rev', ['aTgC', 'AtGc']) # ['CgTa', 'cGtA']
```

### `process_prot(operation, seqs)`
Accepts two variables:
- `operation` defines defines the operation to be applied
- `seqs` list of strings with sequences to process
Returns list of perocessed sequences.

**Valid sequence** should contain 1-letter symbols (case insensetive) of 20 common amino acids ('U' for selenocysteine and 'O' for pyrrolysine doesn't allowed).

**Operations**:
- `lengths` return list with numbers of AA in each sequence(s)
- `molw` return list of protein molecular weight (use the average molecular weight of AA, 110 Da)
- `iso` return list of approximate isoelectric point of given amino acids sequence
- `gravy` return list of GRAVY (grand average of hydropathy) values
- `rewrite` return list of sequences in 3-letter AA code (AA separated by hyphens)

**Examples**
```python
process_prot('iso', ['ACGTWWA', 'ilattwp']) # [5.8, 6.0]
process_prot('rename', ['ACGTwwa']) # ['Ala-Cys-Gly-Thr-Trp-Trp-Ala']
```

### `filter_fastq(input_path, gc_bounds, len_bounds, quality_threshold, output_filename)`
Creates from fastq file new fastq file (in 'fastq_filtrator_resuls' subfolder) with filtered sequences based on specified conditions:
- GC-content, inside interval include borders, or, if single value, not bigger than specified
- length, inside interval include borders, or, if single value, not bigger than specified
- average phred scores, not less than specified

**Options:**
- `input_path` path to the fastq file to be processed (function is case insensetive)
- `gc_bounds` values for GC-content filter.
- `len_bounds` values for lenght filter
- `quality_threshold` value for phred scores filter (`float` or `int`)
- `output_filename` you can specified (defaults is input filename)

**Intervals** for `gc_bounds` and `len_bounds` must be specified (*with bounds are included*):
- as tuple `(lower_bound, upper_bound)`

or

- as single value (`float` or `int`) for upper bound only

**Defaults filters parametrs are**:
- `gc_bounds = (20, 80)`
- `len_bounds = (0, 2**32)`
- `quality_threshold = 0`

**Examples**

```python
"""file.fastq
@A
GTGCATGCGTGCCTGC
+
;AD8@.=BB<7:D.B
@B
ATATTA
+
/++*==
"""

filter_fastq('file.fastq', gc_bounds = 10, output_filename = 'filtered')

"""fastq_filtrator_resuls/filtered.fastq
@B
ATATTA
+
/++*==
"""
```

## `bio_files_processor.py`

It includes tree functions:
- `convert_multiline_fasta_to_oneline` remove line breaks within sequences
- `select_genes_from_gbk_to_fasta` select genes from gbk-file relative to given ones and according to the specified ranges
- `parse_blast_output` select one best hit for each query from file with results of Blast

### `convert_multiline_fasta_to_oneline(input_fasta, output_fasta)`
Convert sequences in fasta files from multiple lines entries with line breaks to single line entries.

Default output filename `one_line_seqs.fasta`.

### `select_genes_from_gbk_to_fasta(input_gbk, genes, n_before, n_after, output_fasta)`
You can pass a list of GOI names as `genes`, specify the up- and downstream ranges around them as `n_before` and `n_after`, and you will get a new fasta-file consists of translations of corresponding regions.
If intervals overlap, regions will be repeated.

Default output filename `results_from_gbk.fasta`.

### `parse_blast_output(input_file, output_file)`
Takes file with results of Blast, select from it one top hit for each queryand save them to txt file sorted alphabetically.

Default output filename `best_Blast_results.txt`.

## Autors
- Kirill Petrikov - main development

Parts of proteins processig code (`prot_seq_tool.py`): Yury Popov, Gulgaz Muradova.    
