import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty


def split_meme_and_cutoff_files(meme_file_path, cutoffs_file_path, motifs_per_file, done_path, argv='no_argv'):

    verify_file_is_not_empty(meme_file_path)
    verify_file_is_not_empty(cutoffs_file_path)

    splitted_meme_dir = os.path.join(os.path.split(meme_file_path)[0], 'memes')
    os.makedirs(splitted_meme_dir, exist_ok=True)

    splitted_cutoffs_dir = os.path.join(os.path.split(cutoffs_file_path)[0], 'cutoffs')
    os.makedirs(splitted_cutoffs_dir, exist_ok=True)

    logger.info(f'{datetime.datetime.now()}: splitting pssms and cuttoffs to:\n'
                f'{splitted_meme_dir}\n'
                f'{splitted_cutoffs_dir}')

    with open(meme_file_path) as meme_f:
        meta_info = ''
        data = ''
        motif_number = 0
        split_number = 0
        add_meta_info = True
        for line in meme_f:
            if add_meta_info:
                if "MOTIF" not in line:
                    meta_info += line
                    continue
                else:
                    add_meta_info = False
            if line.startswith("MOTIF"):
                if motif_number == motifs_per_file:
                    with open(f'{splitted_meme_dir}/{str(split_number).zfill(2)}.txt', 'w') as f:
                        f.write(meta_info + data)
                    data = ''
                    motif_number = 0
                    split_number += 1
                motif_number += 1
            data += line
        # don't forget last batch!!
        with open(f'{splitted_meme_dir}/{str(split_number).zfill(2)}.txt', 'w') as f:
            f.write(meta_info + data)

    with open(cutoffs_file_path) as cutoffs_f:
        data = ''
        motif_number = 0
        split_number = 0
        for line in cutoffs_f:
            if line.startswith("###"):
                if motif_number == motifs_per_file:
                    with open(f'{splitted_cutoffs_dir}/{str(split_number).zfill(2)}.txt', 'w') as f:
                        f.write(data)
                    data = ''
                    motif_number = 0
                    split_number += 1
                motif_number += 1
            data += line
        # don't forget last batch!!
        with open(f'{splitted_cutoffs_dir}/{str(split_number).zfill(2)}.txt', 'w') as f:
            f.write(data)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_file_path', help='A path to a meme file')
    parser.add_argument('cutoffs_file_path', help='A path to a cutoffs file')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--motifs_per_file', type=int, default=5, help='How many motifs will be in each splitted file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging

    if args.verbose:
      logging.basicConfig(level=logging.DEBUG)
    else:
      logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    split_meme_and_cutoff_files(args.meme_file_path, args.cutoffs_file_path, args.motifs_per_file,
                                args.done_file_path, sys.argv)

