from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, load_fasta_to_dict
import logging
logger = logging.getLogger('main')


def remove_sparse_columns(msa_path, out_path, maximal_gap_frequency_allowed_per_column):
    logger.info(f'Removing sparse columns from {msa_path} (allowing columns with gap frequency lower than {maximal_gap_frequency_allowed_per_column})')
    verify_file_is_not_empty(msa_path)

    header_to_sequence, number_of_sequences, msa_length = load_fasta_to_dict(msa_path)
    cleaned_header_to_sequence = dict.fromkeys(header_to_sequence, '')
    for j in range(msa_length):
        column_j = [header_to_sequence[header][j] for header in header_to_sequence]
        gap_frequency = column_j.count('-')/number_of_sequences
        if gap_frequency <= maximal_gap_frequency_allowed_per_column:
            # not a sparse column
            for header in header_to_sequence: # add j'th column
                cleaned_header_to_sequence[header] += header_to_sequence[header][j]
        else:
            logger.debug(f'Removing column #{j}: {column_j}')

    with open(out_path, 'w') as f:
        for header in cleaned_header_to_sequence:
            f.write(f'>{header}\n{cleaned_header_to_sequence[header]}\n')

    logger.info(f'Shortened from {msa_length} to {len(cleaned_header_to_sequence[header])} columns')


if __name__ == '__main__':
    from sys import argv

    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('msa_path', help='A path fasta file with a multiple sequence alignment')
    parser.add_argument('cleaned_msa_path', help='A path to a file in which the filtered MSA will be written to')
    parser.add_argument('--maximal_gap_frequency_allowed_per_column', default=0.1, help='Minimal sequence similarity threshold required',
                        type=lambda x: float(x) if 0 < float(x) < 1
                                                else parser.error(f'The threshold of the maximal gap frequency allowed per column should be between 0 to 1'))
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    remove_sparse_columns(args.msa_path, args.cleaned_msa_path,
                          args.maximal_gap_frequency_allowed_per_column)


