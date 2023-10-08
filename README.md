# DAYWWYD (ver. A)
## Do Anything You Want With Your Data (Almost)
### Expect a new version to be released: L (Literally). Maybe.
> *This is the repo for the homework of the BI Python 2023 course*

DAYWWYD is script for different bioinformatics operations with sequencecs.
It includes three functions:
- `process_na` - operations on nucleic acids sequences
- `process_prot` - operations on proteins (amino acids sequenses)
- `filter_fastq` - filtering sequences from fastq files

## Requirements and dependencies

Each of that functions requires an service module to work, which should be located in the ./src folder.
- `process_na` depend on `na_seq_tool.py`
- `process_prot` depend on `prot_seq_tool.py`
- `filter_fastq` depend on `fastq_filter.py`

Functions work independently of one another and of modules of other functions.

## Functions descriptions
### `process_na(option, seqs)`

Accepts two variables:
- `option`: defines defines the operation to be applied
- `seqs`: list of strings with sequences to process
Returns list of perocessed sequences, the letters case is preserved.

**Valid characters** for sequences: `A, T, G, C, U, a, t, g, c, u`. Chimeric DNA-RNA sequences, i.e. including both `T, t`, and `U, u` aren't allowed.

**Options**:
- `trans` - returns the transcribed sequence. If the original sequence contained `U, u` then they will be replaced by `T, t` as for reverse transcription
- `rev` - returns the inverted sequence, from backward to forward
- `comp` - returns the complementary sequence
- `revcomp` - returns the inverted complementary sequence

**Examples**
```python
process_na('trans', ['ATG']) # [AUG]
process_na('rev', ['aTgC', 'AtGc']) # ['CgTa', 'cGtA']
```

### `process_prot(option, seqs)`
Accepts two variables:
- `option`: defines defines the operation to be applied
- `seqs`: list of strings with sequences to process
Returns list of perocessed sequences.

**Valid sequence** should contain 1-letter symbols (case insensetive) of 20 common amino acids ('U' for selenocysteine and 'O' for pyrrolysine doesn't allowed).

**Options**:
- `lengths`: return list with numbers of AA in each sequence(s)
- `molw`:return list of protein molecular weight (use the average molecular weight of AA, 110 Da)
- `iso`: return list of approximate isoelectric point of given amino acids sequence
- `gravy`: return list of GRAVY (grand average of hydropathy) values
- `rewrite`: return list of sequences in 3-letter AA code (AA separated by hyphens)

**Examples**
```python
process_prot('iso', ['ACGTWWA', 'ilattwp']) # [5.8, 6.0]
process_prot('rename', ['ACGTwwa']) # ['Ala-Cys-Gly-Thr-Trp-Trp-Ala']
```

### `filter_fastq(seqs, gc_bounds, len_bounds, quality_threshold)`
Filters out sequences that satisfy the specified conditions:
- GC-content, inside interval include borders, or, if single value, not bigger than specified
- length, inside interval include borders, or, if single value, not bigger than specified
- average phred scores, not less than specified

**Accepts four variables**:
- `seqs`: dictionary with sequences, names and phred scores
- `gc_bounds`: values for GC-content filter
- `len_bounds`: values for lenght filter
- `quality_threshold`: value for phred scores filter (`float` or `int`)
Returns dictionary with entrys that satisfy the spesified conditions.

**Input fastq dictionary** must be of the form: `['sequence title'] = ('sequense', 'phred scores')`.
Function is case insensetive.

**Intervals** `gc_bounds` and `len_bounds` must be:
- as tuple `(lover_bound, upper_bound)`
- as single value (`float` or `int`) for upper bound

*Intervals bounds are included in filter*

**Defaults filters parametrs are**:
- `gc_bounds = (20, 80)`
- `len_bounds = (0, 2**32)`
- `quality_threshold = 0`

**Examples**

```python
fastq = {"@A": ("GTGCATGCGTGCCTGC", ";AD8@.=BB<7:D.B"),
         "@B":("ATATTA", "/++*==")
        }

filter_fastq(fastq, gc_bounds = 10)
### {"@B":("ATATTA", "/++*==")}
```

## Autors
- Kirill Petrikov - main development

Also parts of proteins processig code (`prot_seq_tool.py`): Yury Popov, Gulgaz Muradova.





    
