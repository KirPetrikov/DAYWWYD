""""Filters out sequences that satisfy the specified conditions:
    - GC-content, inside interval include borders, or, if single value, not bigger than specified
    - length, inside interval include borders, or, if single value, not bigger than specified
    - average phred scores, not less than specified

    Defaults: gc_bounds = (20, 80)
              length_bounds = (0, 2 ** 32)
              quality_threshold = 0

    GATeC v1.0 - that gates opens only for fancy data
"""


def parse_intervals(interval) -> tuple:  # можно добавить проверку неотриц. и min<max
    """Check interval input: int or tuple
       and return values as tuple
    """
    if isinstance(interval, int) or isinstance(interval, float):
        upper_bound = interval
        lover_bound = 0
    else:
        upper_bound = interval[1]
        lover_bound = interval[0]
    bounds = (lover_bound, upper_bound)
    return bounds


def calc_gc_content(seq):
    gc_count = {"g": 0, "G": 0, "c": 0, "C": 0}
    for nuc in seq:
        if nuc in gc_count:
            gc_count[nuc] += 1
    gc_sum = gc_count["g"] + gc_count["G"] + gc_count["c"] + gc_count["C"]
    return 100 * gc_sum / len(seq)


def is_seq_pass_gc_filter(seq_to_check_gc: str, gc_min: float, gc_max: float) -> bool:
    is_lower_pass = calc_gc_content(seq_to_check_gc) >= gc_min
    is_upper_pass = calc_gc_content(seq_to_check_gc) <= gc_max
    return is_lower_pass and is_upper_pass


def is_seq_pass_len_filter(seq_to_check: str, len_min: int, len_max: int) -> bool:
    seq_len = len(seq_to_check)
    is_lower_pass = seq_len >= len_min
    is_upper_pass = seq_len <= len_max
    return is_lower_pass and is_upper_pass


def calc_mean_phred(fastq_phred: str):
    phred_scores = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, '\'': 6, '(': 7, ')': 8, '*': 9, '+': 10,
                    ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
                    '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
                    '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40}
    scores_sum = 0
    for i in fastq_phred:
        scores_sum += phred_scores[i]
    return scores_sum / len(fastq_phred)


def is_seq_pass_phred_filter(fastq_to_check_phred: str, crit_phred: float):
    return calc_mean_phred(fastq_to_check_phred) >= crit_phred
