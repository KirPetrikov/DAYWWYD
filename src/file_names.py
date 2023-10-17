"""Constructs a file name based on
    the requirements provided
"""

from os import path


def file_name_creator(input_path: str, output_file_name: str = '', file_extension: str = '') -> str:
    """Constructs a filename:
        - returns only the name from the full path
        - creates a new name
        - changes extension
    """
    if output_file_name == '':
        if file_extension == '':
            file_name = path.basename(input_path)
        else:
            file_name = '.'.join(path.basename(input_path).split('.')[:-1]) + '.' + file_extension
            # work in case when input file name includes '.'
    else:
        if file_extension == '':
            file_name = output_file_name
        else:
            file_name = output_file_name + '.' + file_extension
    return file_name
