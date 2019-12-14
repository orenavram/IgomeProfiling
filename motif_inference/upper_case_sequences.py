import logging
from Auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, load_fasta_to_dict

def convert_sequences_to_upper(fasta_file, output_path):

    verify_file_is_not_empty(fasta_file)

    header_to_sequence, number_of_sequences, msa_length = load_fasta_to_dict(fasta_file)
    with open(output_path, 'w') as f:
        for header in header_to_sequence:
            f.write(f'>{header}\n{header_to_sequence[header].upper()}\n')

    verify_file_is_not_empty(output_path)


if __name__ == '__main__':
    from sys import argv

    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help='A fasta file')
    parser.add_argument('output_path', help='A fasta file with all the sequences in upper case letters')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    convert_sequences_to_upper(args.fasta_file, args.output_path)
