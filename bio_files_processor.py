"""Provides operations with fasta and gbk files
   Includes two functions:
   - conversion multiline fasta to oneline
   - selection of up- and downstream genes
     from gbk-file relative to given ones
     and according to the specified ranges
"""

from os import path


def file_name_creator(input_path: str = '',
                      output_file_name: str = '',
                      file_extension: str = '') -> str:
    """Constructs a filename:
        - returns only the name from the full path
        - creates a new name
        - changes file extension
    """
    if output_file_name == '':  # use input filename
        if file_extension == '':
            file_name = path.basename(input_path)
        else:
            file_name = '.'.join(path.basename(input_path).split('.')[:-1]) + '.' + file_extension
            # work in case when input file name includes '.'
    else:
        if file_extension == '':
            file_name = output_file_name
        else:
            file_name = output_file_name + '.' + file_extension
    return file_name


def parse_gbk_to_list(path_to_file: str) -> list:
    """Parse gbk-file by CDS in order of appearance
       to list of lists as:
       [
        ['First CDS coordinates',
         'gene (if available) or CDS coordinates',
         'translation (if available) or empty'],
        ['Second CDS', 'gene/CDS', 'translation'],
        etc.
       ]
    """
    with open(path_to_file) as file:
        result_gbk_list = []
        line = file.readline()
        counter = 0
        while not line.startswith('     CDS'):
            line = file.readline()
        result_gbk_list.append([line.strip('CDS \n'), line.strip('CDS \n'), ''])
        for line in file:
            if line.strip().startswith('/gene'):
                result_gbk_list[counter][1] = line.strip(' /gene="\n')
            if line.strip().startswith('/translation'):
                whole_translation = line.strip(' /translation="\n')
                while not line.strip().endswith('"'):
                    line = file.readline()
                    whole_translation += line.strip(' "\n')
                result_gbk_list[counter][2] = whole_translation
            if line.startswith('     CDS'):
                counter += 1
                result_gbk_list.append([line.strip('CDS \n'), line.strip('CDS \n'), ''])
    return result_gbk_list


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


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str = 'shifted_seq'):
    """Takes single entry fasta and rewrite
        the sequence as circular
        starting from the specified index character
    """
    with open(input_fasta) as fasta_file:
        new_name = fasta_file.readline().strip() + '_shifted_to_' + str(shift)
        input_seq = fasta_file.readline().strip()
    if shift < 0:
        shift += len(input_seq)
    shifted_seq = input_seq[shift:len(input_seq) + 1]
    for idx in range(shift):
        shifted_seq += input_seq[idx]
    output_file = file_name_creator('', output_fasta, 'fasta')
    with open(output_file, mode='w') as shifted_fasta:
        shifted_fasta.write(new_name + '\n')
        shifted_fasta.write(shifted_seq)
