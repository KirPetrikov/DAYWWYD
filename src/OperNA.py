"""OperNA v2.0 (ex dna_rna_tools)
This script perform some simple string operations on nucleotide sequences:
"""


def sequence_check(sequences):
    valid_nucleotides = {"a", "t", "g", "c", "u", "A", "T", "G", "C", "U"}
    valid_string_indexes = []
    non_string_indexes = []
    non_valid_nucleotides_indexes = []
    non_valid_ut_mix_indexes = []
    valid_sequences_indexes = []
    for seq in sequences:
        if isinstance(seq, str):
            valid_string_indexes.append(sequences.index(seq))
        else:
            non_string_indexes.append(sequences.index(seq))
    for i in valid_string_indexes:
        seq_set_tmp = set(sequences[i])
        if not (seq_set_tmp <= valid_nucleotides):
            non_valid_nucleotides_indexes.append(i)
        else:
            if (("T" in seq_set_tmp) or ("t" in seq_set_tmp)) and (("U" in seq_set_tmp) or ("u" in seq_set_tmp)):
                non_valid_ut_mix_indexes.append(i)
            else:
                valid_sequences_indexes.append(i)
    return [valid_sequences_indexes, non_string_indexes, non_valid_nucleotides_indexes, non_valid_ut_mix_indexes]


def make_transcribed_seq(sequence):
    if ("T" in sequence) or ("t" in sequence):
        sequence = sequence.replace('t', 'u').replace('T', 'U')
    else:
        sequence = sequence.replace('u', 't').replace('U', 'T')
    return sequence


def make_reversed_seq(sequence):
    return sequence[::-1]


def make_complement_seq(sequence):
    sequence = sequence.replace('c', '%temp%').replace('g', 'c').replace('%temp%', 'g')
    sequence = sequence.replace('C', "%temp%").replace('G', 'C').replace('%temp%', 'G')
    if ("T" in sequence) or ("t" in sequence):
        sequence = sequence.replace('t', '%temp%').replace('a', 't').replace('%temp%', 'a')
        sequence = sequence.replace('T', "%temp%").replace('A', 'T').replace('%temp%', 'A')
    else:
        sequence = sequence.replace('u', '%temp%').replace('a', 'u').replace('%temp%', 'a')
        sequence = sequence.replace('U', "%temp%").replace('A', 'U').replace('%temp%', 'A')
    return sequence


def make_reverse_complement_seq(sequence):
    return make_reversed_seq(make_complement_seq(sequence))


def run_dna_rna_tools(*args):
    if not ({args[-1]} <= {"transcribe", "reverse", "complement", "reverse_complement", "errors"}):
        print("Error! Invalid operation:", args[-1])
        return
    sequences_indexes = sequence_check(args[:-1])
    sequences_indexes_to_proceed = sequences_indexes[0]
    result_sequences = []
    if args[-1] == "transcribe":
        for i in sequences_indexes_to_proceed:
            result_sequences.append(make_transcribed_seq(args[i]))
    elif args[-1] == "reverse":
        for i in sequences_indexes_to_proceed:
            result_sequences.append(make_reversed_seq(args[i]))
    elif args[-1] == "complement":
        for i in sequences_indexes_to_proceed:
            result_sequences.append(make_complement_seq(args[i]))
    elif args[-1] == "reverse_complement":
        for i in sequences_indexes_to_proceed:
            result_sequences.append(make_reverse_complement_seq(args[i]))
    elif args[-1] == "errors":
        for i in sequences_indexes[1]:
            print("Error!", "Argument No", i + 1, "is not string:", args[i])
        for i in sequences_indexes[2]:
            print("Error!", "Not valid nucleotide symbol in Argument No", i + 1, ":", args[i])
        for i in sequences_indexes[3]:
            print("Error!", "DNA/RNA chimera in Argument No", i + 1, ":", args[i])
    if len(result_sequences) == 1:  # Sorry for this crutch, I'm very embarrassed =(
        return result_sequences[0]
    else:
        return result_sequences
