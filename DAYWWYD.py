"""Provides different operations with bioinformatics data
    using three functions:
    - process_na(): nucleic acids sequences conversion to transcript, inverted,
                    complement or revers complement
    - process_prot(): calculation of protein lengths, molecular weights,
                      isoelectric points and GRAVY values, rewrite 1-letter to 3-letter
    - filter_fastq(): filtering sequences from fastq dict by
                      gc-content, length, phred quality

    Requirements (scripts must be in ./src):
    fastq_filter.py
    na_seq_tool.py
    prot_seq_tool.pys
"""

# импортировать сразу нужные функции from ... import
import src.na_seq_tool as nas
import src.prot_seq_tool as ps
import src.fastq_filter as ff
from src.parse_fastq import parse_fastq
from src.parse_fastq import write_fastq


def process_na(operation: str, seqs: list) -> list:
    """Performs operations on list of nucleotide sequences
    Valid options:
    - trans: returns the transcribed sequence (U, u <-> T, t)
    - rev: returns the inverted sequence from backward to forward
    - comp: returns the complementary sequence
    - revcomp: returns the inverted complementary sequence
    """
    valid_options = {
        "trans": nas.transcribe_seq,
        "rev": nas.reverse_seq,
        "comp": nas.complement_seq,
        "revcomp": nas.reverse_complement_seq,
    }
    if not (isinstance(seqs, list)):
        raise ValueError("Invalid data format!")
    for seq in seqs:
        nas.sequence_check(seq)
    if operation in set(valid_options):
        result_sequences = []
        for seq in seqs:
            result_sequences.append(valid_options[operation](seq))
        return result_sequences
    else:
        raise ValueError("Invalid operation!")


def process_prot(operation: str, seqs: list) -> list:
    """Performs operations on amino acids sequences
    Valid options:
    - gravy: calculate GRAVY values
    - iso: calculate isoelectric points
    - molw: calculate molecular weights
    - lengths: calculate sequences lengths
    - rewrite: rewrite 1-letter sequence to 3-letter sequence
    """
    if not (isinstance(seqs, list)):
        raise ValueError("Invalid data format!")
    ps.check_sequences(seqs)
    valid_options = {
        "gravy": ps.calc_gravy,
        "iso": ps.calc_iso_point,
        "molw": ps.calc_protein_mass,
        "lengths": ps.sequence_length,
        "rewrite": ps.transform_to_three_letters
    }
    if operation in valid_options:
        results = []
        for seq in seqs:
            result_tmp = valid_options[operation](seq.upper())
            results.append(result_tmp)
        return results
    else:
        raise ValueError("Invalid operation!")


def filter_fastq(input_path: str,
                 gc_bounds: int | float | tuple = (20, 80),
                 len_bounds: int | float | tuple = (0, 2 ** 32),
                 quality_threshold: int | float = 0,
                 output_filename: str = None,):
    """"Filters out sequences from fastq file by the specified conditions:
        - GC-content, inside interval include borders, or, if single value, not bigger than specified
        - length, inside interval include borders, or, if single value, not bigger than specified
        - average phred scores, not less than specified
        Creates new file with results in 'fastq_filtrator_resuls' folder
    """
    seqs = parse_fastq(input_path)
    gc_lower = ff.parse_intervals(gc_bounds)[0]
    gc_upper = ff.parse_intervals(gc_bounds)[1]
    len_lower = ff.parse_intervals(len_bounds)[0]
    len_upper = ff.parse_intervals(len_bounds)[1]
    selected_seqs = {}
    for seq_name, (seq, comm, phred) in seqs.items():
        gc_condition_check = ff.is_seq_pass_gc_filter(seq, gc_lower, gc_upper)
        len_condition_check = ff.is_seq_pass_len_filter(seq, len_lower, len_upper)
        phred_condition_check = ff.is_seq_pass_phred_filter(phred, quality_threshold)
        if gc_condition_check and len_condition_check and phred_condition_check:
            selected_seqs[seq_name] = seqs[seq_name]
    if output_filename is None:
        output_filename = input_path
    write_fastq(selected_seqs, output_filename)
    return
