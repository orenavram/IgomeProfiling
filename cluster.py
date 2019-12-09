import subprocess
import os
import logging

def cluster_sequences(fasta_file, output_path, threshhold, word_length, throw_sequences_shorter_than):


    # TODO: module load CD hit
    cmd = f'cd-hit -i {fasta_file} ' \
          f'-o {output_path} ' \
          f'-c {threshhold} ' \
          f'-n {word_length} ' \
          f'-l {throw_sequences_shorter_than}'
    subprocess.call(cmd, shell=True)

    # make sure that there are results and the file is not empty
    with open(output_path) as f:
        if len(f.read(10).strip()) == 0:
            # TODO: write error to a global error file
            logger.error(f'CD-hit failed to cluster {output_path}')
            raise RuntimeError(f'CD-hit failed to cluster {output_path}')

    with open(f'{os.path.split(output_path)[0]}/done.txt', 'w'):
        pass


if __name__ == '__main__':
    from sys import argv

    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help='A fasta file to collapse for unique sequences and their counts')
    parser.add_argument('output_prefix', help='A file prefix in which the 2 result files will use as an output path')
    parser.add_argument('--threshhold', type=int, help='Minimal sequence identity required')
    parser.add_argument('--word_length', choices=['2', '3', '4', '5'], default='2',
                        help='A heuristic of CD-hit. Choose of word size:\n5 for similarity thresholds 0.7 ~ 1.0\n4 for similarity thresholds 0.6 ~ 0.7\n3 for similarity thresholds 0.5 ~ 0.6\n2 for similarity thresholds 0.4 ~ 0.5')
    parser.add_argument('--discard', type=int, default=1, help='Include only sequences longer than <$discard> for the analysis. (CD-hit uses only sequences that are longer than 10 amino acids. When the analysis includes shorter sequences, this threshold should be lowered. Thus, it is set here to 1 by default.)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    cluster_sequences(args.fasta_file, args.output_prefix,
                      args.threshhold, args.word_length, args.discard)
