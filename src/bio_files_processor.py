"""Process files
"""


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    if output_fasta is None:
        output_filename = input_fasta
    output_fasta += '.fastq'
    pass


def select_genes_from_gbk_to_fasta(input_gbk: str, genes, n_before, n_after, output_fasta: str = None):
    if output_fasta is None:
        output_fasta = input_gbk
    pass
