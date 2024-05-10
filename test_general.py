import os
import pytest

from general import (DNASequence,
                     RNASequence,
                     AminoAcidSequence,
                     filter_fastq,
                     run_genscan)


@pytest.fixture
def input_data():
    return DNASequence('ATGC')


def test_dnasequence_transcribe(input_data):
    """
    Test DNASequence.transcribe method
    """
    target_values = 'AUGC'
    values_to_check = input_data.transcribe().sequence
    assert target_values == values_to_check


def test_dnasequence_complement(input_data):
    """
    Test DNASequence.complement method
    """
    target_values = 'TACG'
    values_to_check = input_data.complement().sequence
    assert target_values == values_to_check


def test_dnasequence_gc_content(input_data):
    """
    Test DNASequence.gc_content property
    """
    target_values = 50
    values_to_check = input_data.gc_content
    assert target_values == values_to_check


def test_aminoacidsequence_gravy_aa_values():
    """
    Test that hydropathy values of amino acids are correct
    """
    target_values = {'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
                     'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
                     'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
                     'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5
                     }
    values_to_check = AminoAcidSequence.gravy_aa_values
    assert target_values == values_to_check


@pytest.fixture
def input_file():
    file_path = 'test_fastq.fastq'
    content = ('@SRX42\nATGCATGCATGC\n+\nKKKKKKKKKKKK\n'
               '@SRX777\nATGCATGCATG\n+\nKKKKKKKKKK+')
    with open(file_path, mode='w') as file:
        file.write(content)
    yield file_path
    if os.path.exists(file_path):
        os.remove(file_path)


@pytest.fixture
def output_file():
    file_path = 'filtered.fastq'
    yield file_path
    if os.path.exists(file_path):
        os.remove(file_path)


def test_filter_fastq_length(input_file, output_file):
    """
    Test filter_fastq len_thresholds filter
    """
    filter_fastq(input_file, len_thresholds=11)

    with open(output_file) as file:
        file.readline()
        values_to_check = file.readline().strip()
    target_values = 'ATGCATGCATG'
    assert target_values == values_to_check


def test_filter_fastq_quality(input_file, output_file):
    """
    Test filter_fastq quality_threshold filter
    """
    filter_fastq(input_file, quality_threshold=42)

    with open(output_file) as file:
        file.readline()
        values_to_check = file.readline().strip()
    target_values = 'ATGCATGCATGC'
    assert target_values == values_to_check


def test_filter_fastq_gc(input_file, output_file):
    """
    Test filter_fastq gc_thresholds filter
    """
    filter_fastq('test_fastq.fastq', gc_thresholds=49)

    with open(output_file) as file:
        file.readline()
        values_to_check = file.readline().strip()
    target_values = 'ATGCATGCATG'
    assert target_values == values_to_check


def test_run_genscan_incorrect_input():
    """
    Test that a ValueError is raised when the incorrect exon_cutoff are specified
    """
    with pytest.raises(ValueError) as excinfo:
        run_genscan(sequence='ATGC', organism='Ara', exon_cutoff=5)
    assert str(excinfo.value) == ('Incorrect input of "exon_cutoff": 5! '
                                  'Should be: 1.00, 0.50, 0.25, 0.05, 0.02 or 0.01')

