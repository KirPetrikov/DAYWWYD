#!/usr/bin/env python
# coding: utf-8

""""GATeC v1.0
Leaves sequences that satisfy the specified conditions:
- gc-content, inside interval include borders, or, if single value, not bigger than specified
- length, inside interval include borders, or, if single value, not bigger than specified
- average phred scores, not less than specified
"""

PHRED_SCORES = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, '\'': 6, '(': 7, ')': 8, '*': 9, '+': 10,
                ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
                '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
                '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40}


def philter_seqs(seqs: dict,
                 gc_bounds: tuple = (20, 80),
                 length_bounds: tuple = (0, 2 ** 32),
                 quality_threshold: int = 0) -> dict:
    """
    Check does input sequences pass through specified filters
    """
    selected_keys = []
    for key, value in seqs.items():
        if is_seq_pass_phred_filter(value[1], quality_threshold) and is_seq_pass_len_filter(value[0], length_bounds) and is_seq_pass_gc_filter(value[0], gc_bounds):
            selected_keys.append(key)
    selected_seqs = {}
    if len(selected_keys) != 0:
        for key in selected_keys:
            selected_seqs[key] = seqs[key]
    return selected_seqs


def calc_gc_content(seq):
    gc_count = {"g": 0, "G": 0, "c": 0, "C": 0}
    for nuc in seq:
        if nuc in gc_count:
            gc_count[nuc] += 1
    gc_sum = gc_count["g"] + gc_count["G"] + gc_count["c"] + gc_count["C"]
    return 100 * gc_sum / len(seq)


def calc_mean_phred(fastq_phred: str):
    sum_tmp = 0
    for i in fastq_phred:
        sum_tmp += PHRED_SCORES[i]
    return sum_tmp / len(fastq_phred)


def is_seq_pass_gc_filter(seq_to_check_gc: str, crit_gc):
    if isinstance(crit_gc, int) or isinstance(crit_gc, float):
        if calc_gc_content(seq_to_check_gc) <= crit_gc:
            return True
        else:
            return False
    elif calc_gc_content(seq_to_check_gc) >= crit_gc[0] and calc_gc_content(seq_to_check_gc) <= crit_gc[1]:
        return True
    else:
        return False


def is_seq_pass_len_filter(seq_to_check_len: str, crit_len):
    if isinstance(crit_len, int):
        if len(seq_to_check_len) <= crit_len:
            return True
        else:
            return False
    else:
        if len(seq_to_check_len) >= crit_len[0] and len(seq_to_check_len) <= crit_len[1]:
            return True
        else:
            return False


def is_seq_pass_phred_filter(fastq_to_check_phred: str, crit_phred: int):
    if calc_mean_phred(fastq_to_check_phred) >= crit_phred:
        return True
    else:
        return False
