import argparse
import subprocess
import os
from Auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty


def reconstruct_msa(sequences_file_path, output_file_path):
    import subprocess

    # TODO: module load mafft..
    # --auto Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i and FFT-NS-2, according to data size.
    cmd = f'mafft --auto {sequences_file_path} > {output_file_path}'
    logger.info(f'Starting MAFFT. Executed command is:\n{cmd}')
    subprocess.run(cmd, shell=True)

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(output_file_path)

    with open(f'{os.path.split(output_file_path)[0]}/done.txt', 'w'):
        pass


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('sequences_file_path', help='path to a file with unaligned sequences')
    parser.add_argument('output_file_path', help='path to a file in which the aligned sequences will be written')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    reconstruct_msa(args.sequences_file_path, args.output_file_path)





