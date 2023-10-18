"""Provides operations with fasta and gbk files
   Includes two functions:
   - conversion multiline fasta to oneline
   - selection of up- and downstream genes
     from gbk-file relative to given ones
     and according to the specified ranges
"""

from src.parse_gbk import parse_gbk_to_list
from src.file_names import file_name_creator


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = 'one_line_seqs'):
    """Convert sequences in fasta files
       from multiple lines entry with line breaks
       to single line entry
    """
    output_fasta = file_name_creator(output_file_name=output_fasta, file_extension='fasta')
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


def select_genes_from_gbk_to_fasta(input_gbk: str,
                                   genes: list,
                                   n_before: int = 1,
                                   n_after: int = 1,
                                   output_fasta: str = 'results_from_gbk'):
    """Select from gbk-file and write to fasta-file
       genes with their translation
       from specified ranges around genes
       specified by names.
    """
    output_fasta = file_name_creator(output_file_name=output_fasta, file_extension='fasta')
    gbk_list = parse_gbk_to_list(input_gbk)
    selected_genes_idxs = []
    for gene_to_check in genes:
        for idx, gbk_CDS in enumerate(gbk_list):
            if gene_to_check == gbk_CDS[1]:  # разобраться с проверкой доп. символов
                for idx_tmp in range(idx - n_before, idx):
                    selected_genes_idxs.append(idx_tmp)
                for idx_tmp in range(idx + 1, idx + n_after + 1):
                    selected_genes_idxs.append(idx_tmp)
    with open(output_fasta, mode='w') as result_fasta:
        for selected_idx in selected_genes_idxs:
            result_fasta.write('>' + gbk_list[selected_idx][1] + '\n')
            result_fasta.write(gbk_list[selected_idx][2] + '\n')


def parse_blast_output(input_file: str, output_file: str = 'best_Blast_results'):
    """Takes input_file, select from it
        one top hit for each query
        and save them to output_file.txt
        sorted alphabetically
    """
    with open(input_file) as blast_results:
        best_results = []
        line = blast_results.readline()
        while not line.startswith('Query #'):
            line = blast_results.readline()
        for line in blast_results:
            if line.startswith('Description'):
                first_result = blast_results.readline().split('    ')[0]
                best_results.append(first_result.strip('.'))
    best_results.sort()
    output_file = file_name_creator('', output_file, 'txt')
    with open(output_file, mode='w') as best_results_file:
        for result in best_results:
            best_results_file.write(result + '\n')
