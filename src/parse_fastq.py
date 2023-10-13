"""Read and write fastq files.
   When reading parse to dict as:
   'name' = ('sequence', 'comment', 'phred_scores')
   For writing dict of a similar format is required
"""

from os import makedirs
import os.path


def parse_fastq(path_to_file) -> dict:
    """Parse fastq to dict as
    'name' = ('sequence', 'comment', 'phred_scores')
    """
    fastq_dict = {}
    with open(path_to_file) as fastq:
        for line in fastq:
            fastq_name = line.strip()
            fastq_seq = fastq.readline().strip()
            fastq_comment = fastq.readline().strip()
            fastq_phred = fastq.readline().strip()
            fastq_dict[fastq_name] = (fastq_seq, fastq_comment, fastq_phred)
        return fastq_dict


def write_fastq(fastq_dict: dict, output_filename: str):
    """Parse dict to fastq file.
       Create folder 'fastq_filtrator_resuls'
       and create file 'output_filename.fastq' in it
    """
    output_filename += '.fastq'
    output_folder = 'fastq_filtrator_resuls'
    if not os.path.exists(output_folder):
        makedirs(output_folder)
    output_path = os.path.join('fastq_filtrator_resuls', output_filename)
    with open(output_path, mode='w') as file:
        for seq_name, (seq, comm, phred) in fastq_dict.items():
            file.write(seq_name + '\n')
            file.write(seq + '\n')
            file.write(comm + '\n')
            file.write(phred + '\n')
