"""Process files
"""


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """Convert sequences in fasta files
       from multiple lines entry with line breaks
       to single line entry
    """
    if output_fasta is None:
        output_fasta = input_fasta
    output_fasta += '.fastq'
    with open(input_fasta) as fasta_multiline:
        my_list = [fasta_multiline.readline().strip()]
        seq = ''
        for line in fasta_multiline:
            if not line.startswith('>'):
                seq += line.strip()
            else:
                my_list.append(seq)
                my_list.append(line.strip())
                seq = ''
        my_list.append(seq)
    with open(output_fasta, mode='w') as fasta_oneline:
        for line in my_list:
            fasta_oneline.write(line + '\n')


def select_genes_from_gbk_to_fasta(input_gbk: str, genes, n_before, n_after, output_fasta: str = None):
    if output_fasta is None:
        output_fasta = input_gbk
    pass
