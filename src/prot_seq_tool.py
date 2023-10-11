"""Service module for string operations on amino acids sequence
    Main script: DAYWWYD.py
    Calling function: process_prot()

    Provides:
    - sequence check for valid characters
    - calculation of protein lengths, molecular weights, isoelectric points and GRAVY values
    - rewrite 1-letter sequence to 3-letter sequence

    Ancestor: ProtSeqO v2.0 - just open it
"""


def check_sequences(seqs: list):
    """Check amino acid sequence for valid characters"""
    valid_symbols = {'A',
                     'C',
                     'D',
                     'E',
                     'F',
                     'G',
                     'H',
                     'I',
                     'K',
                     'L',
                     'M',
                     'N',
                     'P',
                     'Q',
                     'R',
                     'S',
                     'T',
                     'V',
                     'W',
                     'Y'}
    if not (isinstance(seqs, list)):
        raise ValueError('Enter valid protein sequence')
    for seq in seqs:
        if (not (isinstance(seq, str))) or (not (set(seq.upper()).issubset(valid_symbols))):
            raise ValueError('Enter valid protein sequence')


def calc_gravy(seq: str) -> float:
    """Calculate GRAVY (grand average of hydropathy) value"""
    gravy_aa_values = {'L': 3.8,
                       'K': -3.9,
                       'M': 1.9,
                       'F': 2.8,
                       'P': -1.6,
                       'S': -0.8,
                       'T': -0.7,
                       'W': -0.9,
                       'Y': -1.3,
                       'V': 4.2,
                       'A': 1.8,
                       'R': -4.5,
                       'N': -3.5,
                       'D': -3.5,
                       'C': 2.5,
                       'Q': -3.5,
                       'E': -3.5,
                       'G': -0.4,
                       'H': -3.2,
                       'I': 4.5}
    gravy_aa_sum = 0
    for amino_ac in seq:
        gravy_aa_sum += gravy_aa_values[amino_ac]
    return round(gravy_aa_sum / len(seq), 3)


def calc_total_charge(charged_amino_ac_numbers_list: dict,
                      ph_value: float) -> float:
    """Calculate the approximate total charge for given pH value"""
    n_terminal_charge = 1 / (1 + 10 ** (ph_value - 8.2))
    c_terminal_charge = -1 / (1 + 10 ** (3.65 - ph_value))
    cys_charge = -charged_amino_ac_numbers_list['C'] / (1 + 10 ** (8.18 - ph_value))
    asp_charge = -charged_amino_ac_numbers_list['D'] / (1 + 10 ** (3.9 - ph_value))
    glu_charge = -charged_amino_ac_numbers_list['E'] / (1 + 10 ** (4.07 - ph_value))
    tyr_charge = -charged_amino_ac_numbers_list['Y'] / (1 + 10 ** (10.46 - ph_value))
    his_charge = charged_amino_ac_numbers_list['H'] / (1 + 10 ** (ph_value - 6.04))
    lys_charge = charged_amino_ac_numbers_list['K'] / (1 + 10 ** (ph_value - 10.54))
    arg_charge = charged_amino_ac_numbers_list['R'] / (1 + 10 ** (ph_value - 12.48))
    total_charge = (n_terminal_charge +
                    c_terminal_charge +
                    cys_charge +
                    asp_charge +
                    glu_charge +
                    tyr_charge +
                    his_charge +
                    lys_charge +
                    arg_charge)
    return total_charge


def calc_iso_point(seq: str) -> float:
    """Calculate approximate isoelectric point"""
    charged_amino_ac_numbers = {'C': 0,
                                'D': 0,
                                'E': 0,
                                'Y': 0,
                                'H': 0,
                                'K': 0,
                                'R': 0
                                }
    for amino_ac in seq:
        if amino_ac in charged_amino_ac_numbers:
            charged_amino_ac_numbers[amino_ac] += 1
    total_charge_tmp = 1
    ph_iso_point = -0.1
    while total_charge_tmp > 0:
        ph_iso_point += 0.1
        total_charge_tmp = calc_total_charge(
            charged_amino_ac_numbers,
            ph_iso_point)
    return round(ph_iso_point, 1)


def transform_to_three_letters(seq: str) -> str:
    """
    Transform 1-letter entry to 3-letter
    separated by hyphens.
    """
    amino_acids = {'A': 'Ala',
                   'R': 'Arg',
                   'N': 'Asn',
                   'D': 'Asp',
                   'V': 'Val',
                   'H': 'His',
                   'G': 'Gly',
                   'Q': 'Gln',
                   'E': 'Glu',
                   'I': 'Ile',
                   'L': 'Leu',
                   'K': 'Lys',
                   'M': 'Met',
                   'P': 'Pro',
                   'S': 'Ser',
                   'Y': 'Tyr',
                   'T': 'Thr',
                   'W': 'Trp',
                   'F': 'Phe',
                   'C': 'Cys'}
    new_protein = []
    for aminoacid in seq:
        new_protein.append(amino_acids[aminoacid.upper()])
    return "-".join(new_protein)


def sequence_length(seq: str) -> int:
    """Counts number of aminoacids in given sequence
    """
    return len(seq)


def calc_protein_mass(seq: str) -> int:
    """Calculate molecular weight using the average
    molecular weight of amino acid - 110 Da
    """
    return len(seq) * 110
