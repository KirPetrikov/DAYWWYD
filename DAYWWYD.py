""""Provides different operations with bioinformatics data:

    Requirements (scripts must be in ./src):
    fastq_filter.py
    na_seq_tool.py
    prot_seq_tool.py
"""

import src.na_seq_tool as nas
import src.prot_seq_tool as ps
import src.fastq_filter as ff


def process_na(option: str, seqs: list) -> list:
    """Performs operations on list of nucleotide sequences
    Valid options:
    - trans: returns the transcribed sequence (U, u <-> T, t)
    - rev: returns the inverted sequence from backward to forward
    - comp: returns the complementary sequence
    - revcomp: returns the inverted complementary sequence
    """
    options_dict = {
        "trans": nas.transcribe_seq,
        "rev": nas.reverse_seq,
        "comp": nas.complement_seq,
        "revcomp": nas.reverse_complement_seq,
    }
    if not (isinstance(seqs, list)):
        raise ValueError("Invalid data format!")
    for seq in seqs:
        nas.sequence_check(seq)
    if option in {"trans", "rev", "comp", "revcomp"}:
        result_sequences = []
        for seq in seqs:
            result_sequences.append(options_dict[option](seq))
        return result_sequences
    else:
        raise ValueError("Invalid operation!")


def process_prot(option: str, seqs: list) -> list:
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
    if option in valid_options:
        results = []
        for seq in seqs:
            result_tmp = valid_options[option](seq.upper())
            results.append(result_tmp)
        return results
    else:
        raise ValueError("Invalid operation!")


def filter_fastq(seqs: dict, gc_bounds: tuple = (20, 80), len_bounds: tuple = (0, 2 ** 32),
                 quality_threshold: int = 0) -> dict:
    """"Filters out sequences that satisfy the specified conditions:
        - GC-content, inside interval include borders, or, if single value, not bigger than specified
        - length, inside interval include borders, or, if single value, not bigger than specified
        - average phred scores, not less than specified
    """
    gc_lower = ff.parse_intervals(gc_bounds)[0]
    gc_upper = ff.parse_intervals(gc_bounds)[1]
    len_lower = ff.parse_intervals(len_bounds)[0]
    len_upper = ff.parse_intervals(len_bounds)[1]
    selected_seqs = {}
    for key, value in seqs.items():
        gc_condition_check = ff.is_seq_pass_gc_filter(value[0], gc_lower, gc_upper)
        len_condition_check = ff.is_seq_pass_len_filter(value[0], len_lower, len_upper)
        phred__condition_check = ff.is_seq_pass_phred_filter(value[1], quality_threshold)
        if gc_condition_check and len_condition_check and phred__condition_check:
            selected_seqs[key] = seqs[key]
    return selected_seqs
