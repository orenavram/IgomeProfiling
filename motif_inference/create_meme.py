# More info regarding MEME file format can be found here: http://meme-suite.org/doc/meme-format.html

import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, load_fasta_to_dict, nnk_table
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


def write_pssm(meme_f, letters, msa_name, column_to_letters_frequency_counter, msa_length, number_of_sequences):

    consensus_sequence = ''.join(
        max(column_to_letters_frequency_counter[column], key=column_to_letters_frequency_counter[column].get)
        for column in column_to_letters_frequency_counter)
    meme_f.write(f'MOTIF {consensus_sequence}_{msa_name}\n'
                 f'letter-probability matrix: '
                 f'alength= {len(letters)} '
                 f'w= {msa_length} '
                 f'nsites= {number_of_sequences}\n')
    for column in column_to_letters_frequency_counter:
        # gaps are not counted so the total number of actual participating sequences can
        # be lower than $number_of_sequences
        number_of_participating_sequences = sum(column_to_letters_frequency_counter[column].values())
        column_distribution_str = ' '.join(f'{count/number_of_participating_sequences}'
                                           for count in column_to_letters_frequency_counter[column].values()) + '\n'
        meme_f.write(column_distribution_str)

    meme_f.write('\n\n')


def create_meme_file(msas_path, meme_path, done_path, minimal_number_of_columns_required, argv='no argv'):

    logger.info(f'{datetime.datetime.now()}: generating a new MEME file at {meme_path}')
    letters = sorted(set(letter.upper() for letter in nnk_table.values()))  # don't differentiate between Q and q...

    meme_f = open(meme_path, 'w')
    # write meme file header
    meme_f.write(f'MEME version 4\n\n'
                 f'ALPHABET= {"".join(letters)}\n\n'
                 f'Background letter frequencies\n'
                 f'{get_background_letters_frequency_str(nnk_table)}\n')

    for msa_name in sorted(os.listdir(msas_path)):  # Sorting pssm in meme files by cluster's rank
        # clusterRank_000_uniqueMembers_72_clusterSize_757849.92.faa
        msa_path = os.path.join(msas_path, msa_name)
        logger.info(f'{datetime.datetime.now()}: writing pssm of {msa_path}')
        # make sure that there are results and the msa file is not empty
        verify_file_is_not_empty(msa_path)
        header_to_sequence, number_of_sequences, msa_length = load_fasta_to_dict(msa_path)
        if msa_length < minimal_number_of_columns_required:
            logger.warning(f'{datetime.datetime.now()}: skipping pssm for {msa_path} with only {msa_length} columns '
                           f'(at least {minimal_number_of_columns_required} is required.')
            continue
        column_to_letters_frequency_counter = get_pssm(header_to_sequence, msa_length, letters)
        write_pssm(meme_f, letters, msa_name, column_to_letters_frequency_counter, msa_length, number_of_sequences)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('msas_path',
                        help='A path to a folder with a multiple sequence alignment to be converted to pssms')
    parser.add_argument('meme_path', help='A path to a new/existing MEME file to add the msa PSSM')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--minimal_number_of_columns_required', default=3, type=int,
                        help='MSAs with less than the number of required columns will be skipped')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    create_meme_file(args.msas_path, args.meme_path, args.done_file_path, args.minimal_number_of_columns_required, sys.argv)


