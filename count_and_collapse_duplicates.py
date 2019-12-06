#######################################################################
# Take a FASTA file (where each sequence is in one line) and
# output a FASTA file where each sequence is unique including number of repeats
# Input: (1)  FASTA file (where each sequence is in one line)
#        (2)  ADD to output header the LibType? (YES|NO) [default=NO]
#        (3)  Normalization method (NONE|RpM)
#             [RpM = 1/(Total reads in sample) x 1,000,000, multiply seq-CopyNumber by RpM]
#/groups/pupko/bialik/trash//first_phase_output/TACGCGAT_UH_8_c/TACGCGAT_EXP_5_10-6-18_S1_ALL_head_1000000.UH_8_c.MistakeAllowed.1.AA.fs
# YES
# RpM >> /groups/pupko/bialik/trash/RpM_Factors.txt

import sys

def filter_reads(fasta_file, out_fasta_file, rpm_factors_file):

    sequences_to_counts = {}
    with open(fasta_file) as f:
        for line in f:
            sequence = f.readline().rstrip()
            sequences_to_counts[sequence] = sequences_to_counts.get(sequence, 0) + 1

    sorted_sequences = sorted(sequences_to_counts, key=sequences_to_counts.get, reverse=True)
    total_sequences_count = sum(sequences_to_counts.values())

    factor = 1
    if rpm_factors_file:
        factor = 1_000_000 / total_sequences_count
        with open(rpm_factors_file, 'a') as f:
            f.write(f'{fasta_file}\t{total_sequences_count}\t{factor}')

    with open(out_fasta_file, 'w') as f:
        for i, seq in enumerate(sorted_sequences):
            counts = sequences_to_counts[seq] * factor
            len_seq = len(seq)
            lib = f'{len_seq}'
            if seq.startswith('C') and seq.endswith('C'):
                lib = f'C{len_seq-2}C'
            f.write(f'>Seq_{i+1}_Lib_{lib}_Len_{len_seq}_counts_{counts}\n{seq}\n')


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help='A fasta file to collapse for unique sequences and their counts')
    parser.add_argument('out_fasta_file', help='A fasta file to write the results')
    parser.add_argument('--rpm', help='Normalize counts to "reads per million" (sequence proportion x 1,000,000)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    filter_reads(args.fasta_file, args.out_fasta_file, args.rpm)
