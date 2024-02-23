from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def make_thresholds(threshold: int | float | tuple) -> tuple:
    """Check thresholds input and convert single value to tuple"""
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
    """Filters out sequences from fastq file by the specified conditions:
        - GC-content, inside interval include borders, or, if single value, not bigger than specified
        - length, inside interval include borders, or, if single value, not bigger than specified
        - average phred scores, not less than specified
        Default output file name 'filtered.fastq'
    """
    records = list(SeqIO.parse(input_path, 'fastq'))

    filtered_1_gc_idxs = []
    filtered_2_len_idxs = []
    filtered_3_phred_idxs = []
    filtered_results = []

    min_gc, max_gc = make_thresholds(gc_thresholds)
    min_len, max_len = make_thresholds(len_thresholds)

    for count, record in enumerate(records):
        gc_percent = gc_fraction(record.seq) * 100
        if min_gc <= gc_percent <= max_gc:
            filtered_1_gc_idxs.append(count)

    for idx in filtered_1_gc_idxs:
        if min_len <= len(records[idx].seq) <= max_len:
            filtered_2_len_idxs.append(idx)

    for idx in filtered_2_len_idxs:
        phred_values = records[idx].letter_annotations['phred_quality']
        if sum(phred_values) / len(phred_values) >= quality_threshold:
            filtered_3_phred_idxs.append(idx)

    for idx in filtered_3_phred_idxs:
        filtered_results.append(records[idx])

    with open(output_path, 'w') as file:
        SeqIO.write(filtered_results, file, 'fastq')
