""""Provides different operations with bioinformatics data
"""

import src.OperNA as nas
import src.ProtSeqO as ps

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

def process_prot(option: str, seqs: list) -> list:
    """Performs operations on amino acids sequences:
    - calculate protein lengths, molecular weights, isoelectric points and GRAVY values
    - rewrite 1-letter sequence to 3-letter sequence
    """
    if not (isinstance(seqs, list)):
        raise ValueError("Invalid data format!")
    ps.check_sequences(seqs)
    valid_options = {
        "gravy": ps.calc_gravy,
        "iso": ps.calc_iso_point,
        "rename": ps.transform_to_three_letters,
        "lengths": ps.sequence_length,
        "molw": ps.calc_protein_mass}
    if option in valid_options:
        results = []
        for seq in seqs:
            result_tmp = valid_options[option](seq.upper())
            results.append(result_tmp)
        return results
    else:
        raise ValueError("Invalid operation!")