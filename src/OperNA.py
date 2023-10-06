"""OperNA v2.0 (ex dna_rna_tools)
This script perform some simple string operations on nucleotide sequences:
"""


def sequence_check(sequence: str):
    valid_nucleotides = {"a", "t", "g", "c", "u", "A", "T", "G", "C", "U"}
    seq_set_tmp = set(sequence)
    if not (seq_set_tmp <= valid_nucleotides):
        raise ValueError("Wrong sequence! Invalid character.")
    if (("T" in seq_set_tmp) or ("t" in seq_set_tmp)) and (("U" in seq_set_tmp) or ("u" in seq_set_tmp)):
        raise ValueError("Wrong sequence! DNA/RNA chimera does not allowed.")


def transcribe_seq(sequence: str) -> str:
    transcription_dict = {
        "a": "a", "A": "A",
        "t": "u", "T": "U",
        "u": "t", "U": "T",
        "g": "g", "G": "G",
        "c": "c", "C": "C"
    }
    return "".join([transcription_dict[nuc] for nuc in sequence])


def reverse_seq(sequence: str) -> str:
    return sequence[::-1]


def complement_seq(sequence: str) -> str:
    if ("T" in sequence) or ("t" in sequence):
        complement_dict = {
            "a": "t", "A": "T",
            "t": "a", "T": "A",
            "g": "C", "G": "c",
            "c": "g", "C": "G"
        }
    else:
        complement_dict = {
            "a": "u", "A": "U",
            "u": "a", "U": "A",
            "g": "C", "G": "c",
            "c": "g", "C": "G"
        }
    return "".join([complement_dict[nuc] for nuc in sequence])


def reverse_complement_seq(sequence: str) -> str:
    return reverse_seq(complement_seq(sequence))


def process_na(option: str, seqs: list) -> list:
    """Function process list of nucleotide sequences
    """
    options_dict = {
        "trans": transcribe_seq,
        "rev": reverse_seq,
        "comp": complement_seq,
        "revcomp": reverse_complement_seq,
    }
    if not (isinstance(seqs, list)):
        raise ValueError('Invalid input!')
    for seq in seqs:
        sequence_check(seq)
    if option in {"trans", "rev", "comp", "revcomp"}:
        result_sequences = []
        for seq in seqs:
            result_sequences.append(options_dict[option](seq))
        return result_sequences
    else:
        raise ValueError('Invalid operation!')
