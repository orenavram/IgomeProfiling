import datetime
from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty
import logging
logger = logging.getLogger('main')

def remove_cysteine(fasta_file, out_fasta_file):
    """
    :param fasta_file: a fasta file with sequences
    :param out_fasta_file: a fasta file with the same sequences but flanking Cysteine is removed
    :return:
    """
    logger.info(f'{datetime.datetime.now()}: removing Cysteine loop from {fasta_file}')

    verify_file_is_not_empty(fasta_file)

    f_in = open(fasta_file)
    f_out = open(out_fasta_file, 'w')
    for header in f_in:
        seq = f_in.readline().rstrip()
        if seq.startswith('C') and seq.endswith('C'):
            seq = f'{seq[1:-1]}'  # remove Cys loop
        f_out.write(f'{header}{seq}\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help='A fasta file with aa sequences')
    parser.add_argument('out_fasta_file', help='A fasta file to write the input sequences without the flanking Cysteine')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    remove_cysteine(args.fasta_file, args.out_fasta_file)
