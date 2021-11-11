from os import path
import sys
import logging
import pandas as pd
if path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)
from motif_inference.generate_weblogo import generate_weblogo as create_weblogo

def extract_mapping(file_path: str):
    mapping = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('MOTIF '):
            parts = line.split(' ')[1].split('_')
            mapping[parts[0]] = '_'.join(parts[1:]).rstrip()
    return mapping


def generate_weblogo(motif: str, mapping: dict, cleaned_path: str, output_path: str):
    # Extract meme name
    faa_file = mapping[motif]
    meme_path = path.join(cleaned_path, faa_file)
    weblogo_output_path = path.join(output_path, f'MOTIF_{motif}_{faa_file}')
    logger.info(f'Generating weblogo for {motif}')
    create_weblogo(meme_path, weblogo_output_path, motif)


def extract_weblogos(memes_path, cleaned_path, output_path, done_file_path, motifs_name_path, argv):
    mapping = extract_mapping(memes_path)
    motifs = []
    if motifs_name_path is None:
        motifs = mapping.keys()
    else:
        df_motifs_names = pd.read_csv(motifs_name_path)
        motifs = list(df_motifs_names.columns)[1:] # not take the sample name column.
    
    for motif in motifs:
        generate_weblogo(motif, mapping, cleaned_path, output_path)

    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == "__main__":    
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_file_path', type=str, help='A path to a meme file')
    parser.add_argument('cleaned_path', type=str, help='A csv file with data matrix to model')
    parser.add_argument('output_path', type=str, help='A csv file with data matrix to model')
    parser.add_argument('done_file_path', type=str, help='A path to a file that signals that the module finished running successfully')
    parser.add_argument('--motifs_name_path', type=str, default=None, help='Path for csv file that contain the names of motifs to generate weblogo')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger('main')

    extract_weblogos(args.meme_file_path, args.cleaned_path, args.output_path, args.done_file_path, args.motifs_name_path, sys.argv)
