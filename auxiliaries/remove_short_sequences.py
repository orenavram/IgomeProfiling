import logging
import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, load_fasta_to_dict

def remove_sequences_shorter_than(in_fasta_file, out_fasta_file, minimal_length, argv='no_argv'):

    logger.info(f'{datetime.datetime.now()}: removing all shoter-than-{minimal_length} sequences in {in_fasta_file}')

    verify_file_is_not_empty(in_fasta_file)

    header_to_sequence, number_of_sequences, msa_length = load_fasta_to_dict(in_fasta_file)
    with open(out_fasta_file, 'w') as f:
        for header in header_to_sequence:
            if len(header_to_sequence[header]) >= minimal_length:
                f.write(f'>{header}\n{header_to_sequence[header].upper()}\n')

    verify_file_is_not_empty(out_fasta_file)


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('in_fasta_file', help='A fasta file')
    parser.add_argument('out_fasta_file', help='A fasta file with all the sequences in upper case letters')
    parser.add_argument('min_length', type=int, help='Minimal allowed length')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    remove_sequences_shorter_than(args.in_fasta_file, args.out_fasta_file, args.min_length, sys.argv)
