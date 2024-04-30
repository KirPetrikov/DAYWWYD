"""
The module includes:

- collection of classes fot biological sequences

- filter_fastq : function to filter fastq-files
        by GC-content, length and phred-scores

- telegram_logger : the decorator
        allows you to use a telegram bot to track
        the execution of the decorated function

- run_genscan : API function
        for Genscan Web Server for exons prediction
"""

import datetime
import pandas as pd
import re
import requests
import sys

from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import GC
from dataclasses import dataclass
from dotenv import load_dotenv
from io import BytesIO, StringIO
from os import getenv
from typing import List, TextIO


class InvalidSequenceSymbolError(ValueError):
    """Custom error for BiologicalSequence descendant classes"""
    pass


class BiologicalSequence(ABC):
    """Abstract class to set interface for biological sequences"""

    @abstractmethod
    def check_sequence(self):
        pass

    @abstractmethod
    def __getitem__(self, idx):
        pass

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Class for nucleic acid

    `nucleotides_alphabet` and `complement_dict`
    should be implemented in descendant class

    Methods

    check_sequence : Checks whether characters in a sequence match
        the specified symbols alphabet
    complement : Return complement sequence
    gc_content [property] : Calculates GC-content
    """

    complement_dict = None
    nucleotides_alphabet = None

    def __init__(self, sequence: str):
        self.sequence = sequence

    def complement(self):
        if self.complement_dict is None:
            raise NotImplementedError('It is a basic NA class. '
                                      'You should implement complementarity dictionary '
                                      'in descendant class, e.g. DNASequence.')

        complement_seq = type(self)(''.join([self.complement_dict[nuc] for nuc in self.sequence]))
        return complement_seq

    @property
    def gc_content(self) -> float:
        if not self.check_sequence():
            raise InvalidSequenceSymbolError('Check nucleotide symbols!')

        gc_sum = self.sequence.upper().count('G') + self.sequence.upper().count('C')
        return round(100 * gc_sum / len(self.sequence), 1)

    def check_sequence(self) -> bool:
        if self.nucleotides_alphabet is None:
            raise NotImplementedError('It is a basic NA class. '
                                      'You should implement nucleotides alphabet '
                                      'in descendant class, e.g. DNASequence.')

        seq_symbols = set(self.sequence.upper())
        return seq_symbols.issubset(self.nucleotides_alphabet)

    def __getitem__(self, idx: int):
        return self.sequence.__getitem__(idx)

    def __len__(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return str(self.sequence)

    def __repr__(self) -> str:
        return f'{type(self).__name__}({self.sequence})'


class DNASequence(NucleicAcidSequence):
    """
    Class for DNA sequences

    Methods

    transcribe : Returns `RNASequence` instance with transcribe sequence
    """

    nucleotides_alphabet = {'A', 'C', 'G', 'T'}

    transcription_dict = {'a': 'a', 'A': 'A',
                          't': 'u', 'T': 'U',
                          'g': 'g', 'G': 'G',
                          'c': 'c', 'C': 'C'
                          }

    complement_dict = {'a': 't', 'A': 'T',
                       't': 'a', 'T': 'A',
                       'g': 'c', 'G': 'C',
                       'c': 'g', 'C': 'G'
                       }

    def __init__(self, sequence: str):
        super().__init__(sequence)

    def transcribe(self):
        transcribed_seq = RNASequence(''.join([self.transcription_dict[nuc] for nuc in self.sequence]))
        return transcribed_seq


class RNASequence(NucleicAcidSequence):
    """
    Class for RNA sequences
    """

    nucleotides_alphabet = {'A', 'C', 'G', 'U'}

    complement_dict = {'a': 'u', 'A': 'U',
                       'u': 'a', 'U': 'A',
                       'g': 'c', 'G': 'C',
                       'c': 'g', 'C': 'G'
                       }

    def __init__(self, sequence: str):
        super().__init__(sequence)


class AminoAcidSequence(BiologicalSequence):
    """
    Class for amino acids sequences

    Methods

    check_sequence : Checks whether characters in a sequence match
        the specified symbols alphabet

    gravy [property] : Calculate GRAVY (grand average of hydropathy) value
    """

    aa_alphabet = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
                   }

    gravy_aa_values = {'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
                       'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
                       'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
                       'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5
                       }

    def __init__(self, sequence: str):
        self.sequence = sequence

    def check_sequence(self) -> bool:
        seq_symbols = set(self.sequence.upper())
        return seq_symbols.issubset(self.aa_alphabet)

    @property
    def gravy(self) -> float:
        if not self.check_sequence():
            raise InvalidSequenceSymbolError('Check amino acids symbols!')

        gravy_aa_sum = sum(self.gravy_aa_values[amino_ac] for amino_ac in self.sequence.upper())
        return round(gravy_aa_sum / len(self.sequence), 3)

    def __getitem__(self, idx):
        return self.sequence.__getitem__(idx)

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return str(self.sequence)

    def __repr__(self):
        return f'{type(self).__name__}({self.sequence})'


def make_thresholds(threshold: int | float | tuple) -> tuple:
    """
    Check threshold inputs and convert single value to tuple

    Used in: filter_fastq
    """

    if isinstance(threshold, int) or isinstance(threshold, float):
        lower = 0
        upper = threshold
    else:
        lower = threshold[0]
        upper = threshold[1]
    return lower, upper


def filter_fastq(input_path: str,
                 gc_thresholds: int | float | tuple = (20, 80),
                 len_thresholds: int | float | tuple = (0, 2 ** 32),
                 quality_threshold: int | float = 0,
                 output_path: str = 'filtered.fastq'):
    """
    Filters out sequences from fastq file by the specified conditions:
        - GC-content, inside interval include borders, or, if single value, not bigger than specified
        - length, inside interval include borders, or, if single value, not bigger than specified
        - average phred scores, not less than specified
        Default output file name 'filtered.fastq'

    Params
    ------
    input_path : str
        Path to input fastq-file
    gc_thresholds : int, float or tuple, default (20, 80)
    len_thresholds : int, float or tuple, default (0, 2 ** 32)
    quality_threshold : int or float, default 0
    output_path : str, default 'filtered.fastq'
        Path to output filtered fastq-file
    """

    records_handle = SeqIO.parse(input_path, 'fastq')

    filtered_results = []

    min_gc, max_gc = make_thresholds(gc_thresholds)
    min_len, max_len = make_thresholds(len_thresholds)

    for record in records_handle:
        gc_percent = GC(record.seq)
        phred_values = record.letter_annotations['phred_quality']

        check_gc = min_gc <= gc_percent <= max_gc
        check_len = min_len <= len(record.seq) <= max_len
        check_qual = sum(phred_values) / len(phred_values) >= quality_threshold

        if all((check_gc, check_len, check_qual)):
            filtered_results.append(record)

    with open(output_path, 'w') as file:
        SeqIO.write(filtered_results, file, 'fastq')


def format_time_delta(time_delta: datetime.timedelta) -> str:
    """
    Remove microseconds from the object 'datetime.timedelta'
    if timedelta more than 1 day

    Used in: telegram_logger()
    """

    if time_delta.days > 0:
        time_delta = datetime.timedelta(days=time_delta.days, seconds=time_delta.seconds, microseconds=0)

    return str(time_delta)


def telegram_logger(chat_id: int):
    """
    The decorator allows you to monitor the execution of a function
    and send a message to the telegram bot about completion or an error.
    The message indicates the name of the function, the time spent,
    the error (if any) and attaches a log file.

    A global system variable `TG_API_TOKEN` must be specified
    in `.env` file

    chat_id : int
        user id to whom the message will be sent
    """

    def decorator(func):
        def inner_function(*args, **kwargs):

            load_dotenv()
            tg_api_token = getenv('TG_API_TOKEN')

            try:
                # Streams redirections
                temp_stdout = StringIO()
                temp_stderr = StringIO()
                sys.stdout = temp_stdout
                sys.stderr = temp_stderr

                function_name = func.__name__

                start_time = datetime.datetime.now()
                result = func(*args, **kwargs)

            finally:
                # Save execution time
                end_time = datetime.datetime.now()
                duration = end_time - start_time
                time_text = format_time_delta(duration)

                # Save outputs from stdout and stderr
                log_content = f'Stdout:\n{temp_stdout.getvalue()}\nStderr:\n{temp_stderr.getvalue()}'
                log_content_stream = BytesIO(log_content.encode())

                # Restore system IO-streams
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

                # Save exceptions information
                exc_tuple = sys.exc_info()

                if exc_tuple[0]:  # Error case
                    error_smile = '\U0000274C'
                    message = (f'{error_smile}\n'
                               f'Unfortunately, an error occurred while executing the `{function_name}` after\n'
                               f'`{time_text}` of execution:\n\n'
                               f'`{str(exc_tuple[0])}`\n`{str(exc_tuple[1])}`')

                else:  # Regular case
                    success_smile = '\U00002705'
                    message = (f'{success_smile}\n'
                               f'The `{function_name}` finished successfully in \n`{time_text}`')

                # Set parameters for request
                url = f'https://api.telegram.org/bot{tg_api_token}/sendDocument'
                params = {'chat_id': chat_id,
                          'caption': message,
                          'parse_mode': 'markdown'}
                files = {'document': ('logs.txt', log_content_stream)}

                _ = requests.post(url, params=params, files=files)

            return result
        return inner_function
    return decorator


@dataclass
class GenscanOutput:
    """
    Dataclass containing the results of request processing form Genscan site
    by run_genscan API

    Exons and introns: list of tuples with
        corresponding labels, start and stop coordinates.
    Cds: list of tuples with cds fasta-headers and sequences.
    """

    status: int
    exon_list: list = None
    cds_list: list = None
    intron_list: list = None


def check_input_args(organism: str, exon_cutoff: float, seq_str: str, seq_file: str) -> None:
    """
    Checks input arguments
    Raises error if at least one of them are incorrect

    Used in: run_genscan()
    """

    exon_cutoffs_set = {1.00, 0.50, 0.25, 0.05, 0.02, 0.01}
    organisms_set = {'Vertebrate', 'Arabidopsis', 'Maize'}

    cutoff_check = exon_cutoff not in exon_cutoffs_set
    orgs_check = organism not in organisms_set
    seqs_check = seq_str is None and seq_file is None

    if cutoff_check:
        raise ValueError(f'Incorrect input of "exon_cutoff": {exon_cutoff}! '
                         f'Should be: 1.00, 0.50, 0.25, 0.05, 0.02 or 0.01')
    elif orgs_check:
        raise ValueError(f'Incorrect input of "organism": {organism}! '
                         f'Should be: Vertebrate, Arabidopsis or Maize')
    if seqs_check:
        raise ValueError(f'Incorrect input of sequence: '
                         f'none of sequence string or file name provided! Should specified one of them')


def read_tables(input_data: TextIO) -> List[pd.DataFrame] | None:
    """
    Parse pseudo-table data from plain text

    Returns table with exon labels,
    start and stop coordinates.
    If there are no exons returns None

    Used in: results_parser() -> run_genscan()
    """

    pattern = re.compile(r'\S+')
    exons_tables_list = []

    # Create empty table
    header = ['Exon_ID', 'Type', 'Start', 'End', 'Strand']
    temp_df = pd.DataFrame(columns=header)

    for line in input_data:
        if line == '\n':
            # Checks if two empty line in succession - end of sub-table
            # When functions run - first read line always is data-line,
            # so this condition never True during firs loop
            exons_tables_list.append(temp_df)
            input_data.readline()  # Skip empty line
            line = input_data.readline()
            if line == '\n' or line.startswith('Suboptimal'):
                # Checks if another empty line in succession - end of whole table
                # Or pseudo-table "Suboptimal" begins
                # Both conditions must be checked, since the layout of the results may change
                return exons_tables_list
            else:
                # If no - it's data line for next sub table
                temp_df = pd.DataFrame(columns=header)

        curr_row = re.findall(pattern, line.strip())

        # Checks if there is at least one exon
        if curr_row[0] == 'NO':
            return None

        data_to_add = [curr_row[0], curr_row[1], int(curr_row[3]), int(curr_row[4]), curr_row[2]]
        temp_df.loc[len(temp_df)] = data_to_add

        input_data.readline()  # Skip empty line


def read_cds_sequence(input_data: TextIO, line: str) -> tuple:
    """
    Read part of output with fasta-format sequences
    Returns list of tuples with fasta-header and sequence

    Used in: results_parser() -> run_genscan()
    """

    curr_seq = ''
    seq_header = line.strip()

    for _ in input_data:
        next_line = input_data.readline()
        if next_line == '\n':
            return seq_header, curr_seq
        else:
            curr_seq += next_line.strip()


def results_parser(text_data: str) -> tuple | None:
    """
    Parse results of requests.
    Reads and returns a table with exons, a list of coding sequences

    Used in: run_genscan()
    """
    cds_list = []

    with StringIO(text_data) as results:
        exons_tables = []
        table_start_flag = False
        for line in results:
            # Checks pseudo-table header
            if line.startswith('Gn.Ex'):
                table_start_flag = True
                counter = 0
                while counter < 5:
                    results.readline()
                    counter += 1

            if table_start_flag:
                exons_tables = read_tables(results)
                table_start_flag = False

            # Checks fasta-format header for predicted cds, not peptide
            if line.startswith('>') and 'GENSCAN_predicted_CDS' in line:
                cds_list.append(read_cds_sequence(results, line))

    if len(cds_list) == 0:
        cds_list = None

    return exons_tables, cds_list


def get_introns(exon_data: pd.DataFrame) -> list | None:
    """
    Calculates intron coordinates
    from exons table

    Used in: run_genscan()
    """

    exon_types = {'Init', 'Term', 'Intr'}

    if exon_data['Type'].isin(exon_types).sum() < 2:
        # Checks if there is a single exon and no introns
        return

    ex_start = exon_data.loc[exon_data['Type'].isin(exon_types)].loc[:, 'Start'].argmin()
    ex_stop = exon_data.loc[exon_data['Type'].isin(exon_types)].loc[:, 'Start'].argmax()

    introns = []

    # Strandness detection
    if exon_data.iloc[ex_start, 4] == '+':
        for i in range(ex_start, ex_stop):
            introns.append((exon_data.iloc[i]['Exon_ID'], exon_data.iloc[i, 3] + 1, exon_data.iloc[i + 1, 2] - 1))
    else:
        for i in range(ex_stop, ex_start, -1):
            introns.append((exon_data.iloc[i]['Exon_ID'], exon_data.iloc[i, 3] - 1, exon_data.iloc[i - 1, 2] + 1))

    return introns


def df_to_list(df: pd.DataFrame) -> tuple:
    """
    Convert df with exons to tuple

    Used in: run_genscan()
    """

    return tuple(df.loc[:, ['Exon_ID', 'Start', 'End']].apply(lambda x: tuple(x.tolist()), axis=1).tolist())


def run_genscan(sequence: str = None,
                sequence_file=None,
                organism: str = 'Vertebrate',
                exon_cutoff: float = 1.00,
                sequence_name: str = '') -> GenscanOutput:
    """
    Sends a request to the Genscan Web Server with the specified parameters.
    From the results obtained, it saves a table with exons and the found CDS
    into an GenscanOutput class object.
    Based on the results of the search for exons,
    calculates the coordinates of introns and
    also saves them into the same GenscanOutput instance.

    Params
    ------
    sequence : str
        nucleotide sequence
    sequence_file : str
        path to fasta-file
        For case both text sequence and file are specified
        file are ignored (with warning)
    organism : {'Vertebrate', 'Arabidopsis', 'Maize'}
        default: 'Vertebrate'
    exon_cutoff : {1.00, 0.50, 0.25, 0.05, 0.02, 0.01},
        default: 1.00
    sequence_name : str [optional]

    return : GenscanOutput class object
    """

    check_input_args(organism, exon_cutoff, sequence, sequence_file)

    url = 'http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi'
    hd = {'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:123.0) Gecko/20100101 Firefox/123.0',
          'Referer': 'http://hollywood.mit.edu/GENSCAN.html'
          }

    options_print = 'Predicted CDS and peptides'

    if sequence is None:
        # Send sequence from file
        with open(sequence_file, 'rb') as f:
            data = f.read()
            request_data = {'-o': organism,
                            '-e': exon_cutoff,
                            '-n': sequence_name,
                            '-p': options_print,
                            '-u': data
                            }
    else:
        # Send sequence from text-form
        request_data = {'-o': organism,
                        '-e': exon_cutoff,
                        '-n': sequence_name,
                        '-p': options_print,
                        '-s': sequence
                        }

        if sequence_file is not None:
            print('Warning: both a sequence string and a file are provided. File will be ignored.',
                  file=sys.stderr
                  )

    response = requests.post(url,
                             headers=hd,
                             data=request_data
                             )

    if response.status_code != 200:
        print(f'Warning: response status code is {response.status_code}! Something went wrong.',
              file=sys.stderr
              )

    gs_result = GenscanOutput(response.status_code)

    exons_tables, cds = results_parser(response.text)

    if exons_tables is not None:
        exons = []
        introns = []
        for ex_table in exons_tables:
            introns.extend(get_introns(ex_table))
            exons.extend(df_to_list(ex_table))
            gs_result.exon_list = exons
        if introns[0] is not None:
            gs_result.intron_list = introns

    gs_result.cds_list = cds

    return gs_result

