"""Parse fastq files to dict as:
   'name' = ('sequence', 'comment', 'phred_scores')
"""


def parse_fastq(path_to_file) -> dict:
    fastq_dict = {}
    with open(path_to_file) as fastq:
        for line in fastq:
            fastq_name = line.strip()
            fastq_seq = fastq.readline().strip()
            fastq_comment = fastq.readline().strip()
            fastq_phred = fastq.readline().strip()
            fastq_dict[fastq_name] = (fastq_seq, fastq_comment, fastq_phred)
        return fastq_dict
