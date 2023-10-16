"""Parse gbk-files by CDS to list of lists
   with CDS coordinate,
   gene name (if available),
   translation (if available)
"""


def parse_gbk_to_list(path_to_file: str) -> list:
    """Parse gbk-file by CDS in order of appearance
       to list of lists as:
       [
            ['First CDS coordinates',
            'gene (if available) or CDS coordinates',
            'translation (if available) or empty
            ],
            ['Second CDS', 'gene/CDS', 'translation'],
            etc.
       ]
    """
    with open(path_to_file) as file:
        result_gbk_list = []
        line = file.readline()
        counter = 0
        while not line.startswith('     CDS'):
            line = file.readline()
        result_gbk_list.append([line.strip('CDS \n'), line.strip('CDS \n'), ''])
        for line in file:
            if line.strip().startswith('/gene'):
                result_gbk_list[counter][1] = line.strip(' /gene="\n')
            if line.strip().startswith('/translation'):
                whole_translation = line.strip(' /translation="\n')
                while not line.strip().endswith('"'):
                    line = file.readline()
                    whole_translation += line.strip(' "\n')
                result_gbk_list[counter][2] = whole_translation
            if line.startswith('     CDS'):
                counter += 1
                result_gbk_list.append([line.strip('CDS \n'), line.strip('CDS \n'), ''])
    return result_gbk_list
