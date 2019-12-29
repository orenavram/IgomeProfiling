from weblogolib import *
from corebio.seq import unambiguous_protein_alphabet
import sys
import logging
logger = logging.getLogger('main')


def generate_weblogo(in_file, out_prefix, title):
    logger.info(f'Generating logo for {out_prefix}')
    with open(in_file) as f:
        seqs = read_seq_data(f, alphabet=unambiguous_protein_alphabet)
    data = LogoData.from_seqs(seqs)
    options = LogoOptions()
    options.title = title
    options.resolution = 1500
    options.unit_name = 'probability'
    options.fineprint = ''  # remove "version number bla bla" from the figure
    format = LogoFormat(data, options)
    with open(f'{out_prefix}.png', 'wb') as f:
        f.write(png_formatter(data, format))


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('msa_path', help='A path to a meme file with motifs against which a set of random peptides will be scanned')
    parser.add_argument('output_prefix', help='A prefix (without file extension) to a path in which the weblogo saved')
    parser.add_argument('--title', default="NO TITLE", help='A path to a file that signals that the script was finished running successfully.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    generate_weblogo(args.msa_path, args.output_prefix, args.title)

