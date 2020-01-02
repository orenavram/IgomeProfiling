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

def convert_sequences_to_upper(in_fasta_file, out_fasta_file, done_file_path, argv='no argv'):

    logger.info(f'{datetime.datetime.now()}: upper casing all sequences in {in_fasta_file}')

    verify_file_is_not_empty(in_fasta_file)

    header_to_sequence, number_of_sequences, msa_length = load_fasta_to_dict(in_fasta_file)
    with open(out_fasta_file, 'w') as f:
        for header in header_to_sequence:
            f.write(f'>{header}\n{header_to_sequence[header].upper()}\n')

    verify_file_is_not_empty(out_fasta_file)

    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('in_fasta_file', help='A fasta file')
    parser.add_argument('out_fasta_file', help='A fasta file with all the sequences in upper case letters')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    convert_sequences_to_upper(args.in_fasta_file, args.out_fasta_file, args.done_file_path, sys.argv)
