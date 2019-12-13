# TODO: upper case all the msa!!
# More info regarding MEME file format can be found here: http://meme-suite.org/doc/meme-format.html

import os
import sys
sys.path.insert(0, '/Users/Oren/Dropbox/Projects/gershoni/src/')
print(sys.version)
from Auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, load_msa, nnk_table
import logging
logger = logging.getLogger('main')


def get_pssm(header_to_sequence, msa_length, letters):
    # a dictionary that maps every column to a dictionary that holds a frequency counter of
    # the column letters, i.e., the letters distribution in this column
    column_to_letters_frequency_counter = {column_number: dict.fromkeys(letters, 0)
                                             for column_number in range(msa_length)}

    for header in header_to_sequence:
        sequence = header_to_sequence[header]
        for column_number in range(msa_length):
            current_letter = sequence[column_number]
            if current_letter != '-':  # ignore gaps
                column_to_letters_frequency_counter[column_number][current_letter] += 1

    return column_to_letters_frequency_counter


def get_background_letters_frequency_str(nnk_table):
    aa_to_frequency = {}
    for codon in nnk_table:
        aa = nnk_table[codon].upper()  # don't differentiate between Q and q...
        aa_to_frequency[aa] = aa_to_frequency.get(aa, 0) + 1

    total_number_of_codons = len(nnk_table)
    result = ''
    for aa in sorted(aa_to_frequency):  # go over the amino acids by ABC
        aa_probability = aa_to_frequency[aa]/total_number_of_codons
        result += f'{aa} {aa_probability} '

    return result.rstrip()  # remove last redundant space


def add_pssm_to_meme_file(msa_path, meme_path, add_header):
    if add_header:
        logger.info(f'Generating a new MEME file at {meme_path}')

    logger.info(f'Calculating PSSM of {msa_path}')
    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(msa_path)

    header_to_sequence, number_of_sequences, msa_length = load_msa(msa_path)
    letters = sorted(set(letter.upper() for letter in nnk_table.values()))  # don't differentiate between Q and q...
    column_to_letters_frequency_counter = get_pssm(header_to_sequence, msa_length, letters)

    consensus_sequence = ''.join(max(column_to_letters_frequency_counter[column], key=column_to_letters_frequency_counter[column].get)
                                 for column in column_to_letters_frequency_counter)

    mode = 'a'  # append to an existing file
    meta_info = ''
    if add_header:
        # override previous file!!
        mode = 'w'
        meta_info = f'MEME version 4\n\n' \
                    f'ALPHABET= {"".join(letters)}\n\n' \
                    f'Background letter frequencies\n' \
                    f'{get_background_letters_frequency_str(nnk_table)}\n'
    else:
        # the file already exists and contains at least one PSSM
        # just add some new lines before the next PSSM
        meta_info += '\n\n'
        assert os.path.exists(meme_path), \
            f"add_header parameter wasn't set but as if meme_path exists but it does not!\n{meme_path}\n"

    msa_name = os.path.split(os.path.splitext(msa_path)[0])[1]
    meta_info += f'MOTIF {consensus_sequence}_{msa_name}\n'
    meta_info += f'letter-probability matrix: ' \
                 f'alength= {len(letters)} ' \
                 f'w= {msa_length} ' \
                 f'nsites= {number_of_sequences}\n'

    with open(meme_path, mode) as f:
        f.write(meta_info)
        for column in column_to_letters_frequency_counter:
            # gaps are not counted so the total number of actual participating sequences can
            # be lower than $number_of_sequences
            number_of_participating_sequences = sum(column_to_letters_frequency_counter[column].values())
            column_distribution_str = ' '.join(f'{count/number_of_participating_sequences}'
                                               for count in column_to_letters_frequency_counter[column].values()) + '\n'
            f.write(column_distribution_str)


if __name__ == '__main__':
    from sys import argv

    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('msa_path', help='A path fasta file with a multiple sequence alignment')
    parser.add_argument('meme_path', help='A path to a new/exisiting MEME file to add the msa PSSM')
    parser.add_argument('--header', action='store_true', help='Set to TRUE only when the first PSSM is written (will add a file header according to the MEME format. More info regarding MEME file format can be found here: http://meme-suite.org/doc/meme-format.html)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    add_pssm_to_meme_file(args.msa_path, args.meme_path,
                          True if args.header else False)


