"""
Provides operations with fasta and gbk files

Includes:
- data-class `FastaRecord`
- context manager `OpenFasta`
- functions:
    `convert_multiline_fasta_to_oneline`
    `parse_blast_output`
    `change_fasta_start_pos`
    `select_genes_from_gbk_to_fasta`
"""

from dataclasses import dataclass


@dataclass
class FastaRecord:
    """
    Store fasta-records

    Params
    ------
    id : str
        Sequence accession ID
    seq : str
         Nucleotides/amino acids sequence
    description : str
         Header part after id
    wish_beauty : bool, default True
         Do you want to add beauty?
    """

    id: str
    seq: str
    description: str
    wish_beauty: bool = True
    _beauty: str = ('\n             __\n'
                    '        _   /  |''\n'
                    '       | \  \/_/\n'
                    '       \_\| / __              \n'
                    '          \/_/__\           .--=/~\\\n'
                    '   ____,__/__,_____,______)/  /{~}}}\n'
                    '   -,-----,--\--,-----,---,\  \{{{~}\n'
                    '           __/\_            --=.\}/\n'
                    '          /_/ |\\\\\n'
                    '               \/')

    def __repr__(self):
        if len(self.seq) <= 5:
            print_seq = self.seq
        else:
            print_seq = self.seq[:5] + '...'
        if self.wish_beauty:
            return f"<class FastaRecord>, id='{self.id}', seq='{print_seq}'" + self._beauty
        else:
            return f"<class FastaRecord>, id='{self.id}', seq='{print_seq}'. =("


class FastaFormatError(ValueError):
    """
    Custom error raised, if fasta-file does not start with '>'
    """

    pass


class OpenFasta:
    """
    The context manager to read fasta-file and
    parse them to return FastaRecord class objects.

    Params
    ------
    fasta_path : str
        Path to fasta-file
    wish_beauty : bool, default False
        If you wish!
    """

    def __init__(self, fasta_path: str, wish_beauty: bool = False):
        self.fasta_path = fasta_path
        self.fasta_file = open(fasta_path)
        self.header = self.fasta_file.readline()
        self.wish_beauty = wish_beauty

    def __enter__(self):
        if not self.header.startswith('>'):
            raise FastaFormatError('Invalid fasta-file format: must begin with ">"')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.fasta_file.close()

    def __next__(self):
        if self.header == '':
            raise StopIteration
        header = self.header.strip()[1:]
        fasta_id = header.partition(' ')[0]
        description = header.partition(' ')[2]
        current_seq = ''
        current_line = self.fasta_file.readline()
        while not current_line.startswith('>'):
            current_seq += current_line.strip()
            current_line = self.fasta_file.readline()
            if current_line == '':
                record = FastaRecord(fasta_id, current_seq, description, self.wish_beauty)
                break
        self.header = current_line
        record = FastaRecord(fasta_id, current_seq, description, self.wish_beauty)
        return record

    def __iter__(self):
        return self

    def read_record(self):
        if self.header == '':
            return FastaRecord('', '', '', False)
        return self.__next__()

    def read_records(self):
        full_fasta = []
        for record in self.__iter__():
            full_fasta.append(record)
        return full_fasta

    def __repr__(self):
        return f'{type(self)}. {self.fasta_path}'


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = 'one_line_seqs.fasta'):
    """
    Convert sequences in fasta files
    from multiple lines entry with line breaks
    to single line entry

    Params
    ------
    input_fasta : str
        Path to input fasta-file
    output_fasta : bool, default 'one_line_seqs.fasta'
        Path to output fasta-file
    """

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


def parse_blast_output(input_file: str, output_file: str = 'best_Blast_results.txt'):
    """
    Takes input_file, select from it
    one top hit for each query
    and save them to output_file.txt
    sorted alphabetically

    Params
    ------
    input_file : str
        Path to input txt-file
    output_file : bool, default 'best_Blast_results.txt'
        Path to output txt-file
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

    with open(output_file, mode='w') as best_results_file:
        for result in best_results:
            best_results_file.write(result + '\n')


def change_fasta_start_pos(input_fasta: str, shift_idx: int, output_fasta: str = 'shifted.fasta'):
    """
    Takes single entry fasta and rewrite
    the sequence as circular
    starting from the specified index character.

    Nucleotide position indexing starts from 1,
    means shift = 0 and shift = 1 change nothing

    Params
    ------
    input_fasta : str
        Path to input fasta-file
    shift_idx : int
        Index of nucleotide from which new sequence will start
    output_fasta : str, default 'shifted.fasta'
         Path to output fasta-file
    """
    
    with open(input_fasta) as fasta_file:
        new_name = fasta_file.readline().strip() + '_shifted_to_' + str(shift_idx)
        input_seq = fasta_file.readline().strip()

    if shift_idx < 0:
        shift_idx += len(input_seq)
    elif shift_idx > 1:
        shift_idx -= 1
    else:
        shift_idx = 0

    shifted_seq = input_seq[shift_idx:len(input_seq) + 1]
    shifted_seq += input_seq[:shift_idx]

    with open(output_fasta, mode='w') as shifted_fasta:
        shifted_fasta.write(new_name + '\n')
        shifted_fasta.write(shifted_seq)


def parse_gbk_to_list(path_to_file: str) -> list:
    """
    Parse gbk-file by CDS in order of appearance
       to list of lists as:
       [
        ['First CDS coordinates',
         'gene name (if available) or CDS coordinates',
         'translation (if available) or empty'],
        ['Second CDS coordinates', 'gene/CDS', 'translation'],
        etc.
       ]

    Used in: `select_genes_from_gbk_to_fasta`
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


def select_genes_from_gbk_to_fasta(input_gbk: str,
                                   genes: list | tuple,
                                   n_before: int = 1,
                                   n_after: int = 1,
                                   output_fasta: str = 'results_from_gbk.fasta'):
    """
    Select from gbk-file genes names with their translation
    from specified ranges around genes specified by names.
    Write selected entries to fasta-file

    Params
    ------
    input_gbk : str
        Path to input fasta-file
    genes : list | tuple
        List of genes names for which
        adjacent genes need to be extracted
    n_before: int, default 1
        Range of upstream genes to be extracted
    n_after: int, default 1
        Range of downstream genes to be extracted
    output_fasta : str, default 'results_from_gbk.fasta'
         Path to output fasta-file
    """

    gbk_list = parse_gbk_to_list(input_gbk)

    selected_genes_idxs = []

    for gene_to_check in genes:
        for idx, gbk_CDS in enumerate(gbk_list):
            if gene_to_check == gbk_CDS[1]:
                for idx_tmp in range(idx - n_before, idx):
                    selected_genes_idxs.append(idx_tmp)
                for idx_tmp in range(idx + 1, idx + n_after + 1):
                    selected_genes_idxs.append(idx_tmp)

    with open(output_fasta, mode='w') as result_fasta:
        for selected_idx in selected_genes_idxs:
            result_fasta.write('>' + gbk_list[selected_idx][1] + '\n')
            result_fasta.write(gbk_list[selected_idx][2] + '\n')
