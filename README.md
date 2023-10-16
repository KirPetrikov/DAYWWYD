# DAYWWYD (ver. B)
## Do Anything You Want With Your Data (Better than ver. A)
*Expect a new version to be released: L (Literally). Maybe.*
> *This is the repo for the homework of the BI Python 2023 course*

**DAYWWYD is set of scripts for different bioinformatics operations with sequences.**

It consit of two main scrits:
- `DAYWWYD.A`
- `bio_files_processor.py`

## Requirements and dependencies
`DAYWWYD.A` requires:
- `na_seq_tool.py`
- `prot_seq_tool.py`
- `fastq_filter.py`

`bio_files_processor` requires:
- `parse_gbk.py`

Service modules should be located in the `./src` folder.

## DAYWWYD.A

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
- `output_filename` you can specified (defaults is input file)

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

filter_fastq(fastq, gc_bounds = 10)

"""fastq_filtrator_resuls/filtered.fastq
@B
ATATTA
+
/++*==
"""
```

## bio_files_processor

It includes two functions:
- `convert_multiline_fasta_to_oneline` remove line breaks within sequences
- `select_genes_from_gbk_to_fasta` select genes from gbk-file relative to given ones and according to the specified ranges

### `convert_multiline_fasta_to_oneline(input_fasta, output_fasta)`
Convert sequences in fasta files from multiple lines entry with line breaks to single line entry. Default output file name `one_line_fasta.fasta`.

### `select_genes_from_gbk_to_fasta(input_gbk, genes, n_before, n_after, output_fasta)`
You can pass a list of GOI names as `genes`, specify the up- and downstream ranges around them as `n_before` and `n_after`, and you will get a new fasta-file consists of translations of corresponding regions.

Select from gbk-file and write to fasta-file
       genes with their translation
       from specified ranges around genes
       specified by names.

## Autors
- Kirill Petrikov - main development

Also parts of proteins processig code (`prot_seq_tool.py`): Yury Popov, Gulgaz Muradova.    
