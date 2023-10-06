""""Provides different operations with bioinformatics data
"""

import src.OperNA as nas


def process_na(option: str, seqs: list) -> list:
    """Function process list of nucleotide sequences
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
