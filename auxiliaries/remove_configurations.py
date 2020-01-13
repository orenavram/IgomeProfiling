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

def remove_configurations(in_fasta_file, out_fasta_file, allowed_configurations, argv='no_argv'):

    logger.info(f'{datetime.datetime.now()}: removing all configurations that are not one of these:\n'
                f'{allowed_configurations}\n'
                f'From {in_fasta_file}')

    verify_file_is_not_empty(in_fasta_file)

    header_to_sequence, number_of_sequences, msa_length = load_fasta_to_dict(in_fasta_file)
    with open(out_fasta_file, 'w') as f:
        for header in header_to_sequence:
            for conf in allowed_configurations:
                if f'lib_{conf}_' in header or f'Type_{conf}' in header:
                    f.write(f'>{header}\n{header_to_sequence[header].upper()}\n')
                    break

    verify_file_is_not_empty(out_fasta_file)


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('in_fasta_file', help='A fasta file')
    parser.add_argument('out_fasta_file', help='A fasta file with all the sequences in upper case letters')
    parser.add_argument('allowed_configurations', type=lambda x: x.upper().split(','), help='Which configuration types should be kept')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    remove_configurations(args.in_fasta_file, args.out_fasta_file, args.allowed_configurations, sys.argv)
