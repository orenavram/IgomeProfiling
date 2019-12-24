import datetime
import sys
import logging
logger = logging.getLogger('main')

def count_and_collapse(fasta_file, out_fasta_file, rpm_factors_file, done_file_path, argv='no argv'):
    """
    :param fasta_file: a fasta file with non unique sequences
    :param out_fasta_file: a fasta file with the sequences from the input file but with summarized counts (of
    duplicates), normalized to rpm (if $rpm_factors_file is provided) and collapsed (i.e., without duplicates)
    unique sequences.
    :param rpm_factors_file: a path to which the reads per million (RPM) normalization factor should be written to
    :return:
    """
    open_function = open
    mode = 'r'
    if fasta_file.endswith('gz'):
        import gzip
        open_function = gzip.open
        mode = 'rt'

    logger.info(f'{datetime.datetime.now()}: counting and collapsing {fasta_file} (mode: {mode})')

    sequences_to_counts = get_sequences_frequency_counter(fasta_file, open_function, mode)

    sorted_sequences = sorted(sequences_to_counts, key=sequences_to_counts.get, reverse=True)
    total_sequences_count = sum(sequences_to_counts.values())

    factor = 1
    if rpm_factors_file:
        factor = 1_000_000 / total_sequences_count
        with open(rpm_factors_file, 'a') as f:
            f.write(f'{fasta_file}\t{total_sequences_count}\t{factor}\n')

    with open(out_fasta_file, 'w') as f:
        for i, seq in enumerate(sorted_sequences):
            counts = sequences_to_counts[seq] * factor
            len_seq = len(seq)
            lib = f'{len_seq}'
            if seq.startswith('C') and seq.endswith('C'):
                lib = f'C{len_seq-2}C'
            f.write(f'>seq_{i+1}_lib_{lib}_len_{len_seq}_counts_{counts}\n{seq}\n')

    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')



def get_sequences_frequency_counter(fasta_file, open_function, mode):
    sequences_to_counts = {}
    with open_function(fasta_file, mode) as f:
        for line in f:
            # skip headers. Read only sequences.
            sequence = f.readline().rstrip()
            sequences_to_counts[sequence] = sequences_to_counts.get(sequence, 0) + 1
    return sequences_to_counts


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help='A fasta file to collapse for unique sequences and their counts')
    parser.add_argument('out_fasta_file', help='A fasta file to write the results')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script was finished running successfully.')
    parser.add_argument('--rpm', help='Normalize counts to "reads per million" (sequence proportion x 1,000,000)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    count_and_collapse(args.fasta_file, args.out_fasta_file,
                       args.rpm, args.done_file_path, sys.argv)
