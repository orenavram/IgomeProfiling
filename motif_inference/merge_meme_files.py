import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

import logging
logger = logging.getLogger('main')


#TODO: consider writing pssm sorted accross (biological condition) by their cluster size (currently it's sorted within sample)

def merge_meme_files(motif_inference_path, biological_condition, merged_meme_path, done_path, samples_to_skip, argv='no argv'):
    """
    :param motif_inference_path: A path in which each folder corresponds to a sample and contains a meme file for the
    motifs in this sample.
    :param biological_condition: A biological condition to merge its samples meme files
    :param meme_file_name: An optional file name for the output file
    :return: A merged meme file at $motif_inference_path/$biological_condition/$meme_file_name
    """
    logger.info(f'{datetime.datetime.now()}: merging meme files of {biological_condition}')

    merged_meme_f = open(merged_meme_path, 'w')
    first_meme = True
    for sample_name in sorted(os.listdir(motif_inference_path)):  # sorted by sample name
        dir_path = os.path.join(motif_inference_path, sample_name)
        if not os.path.isdir(dir_path) or biological_condition not in sample_name:
            # skip file or folders of non-related biological condition
            continue
        if sample_name in samples_to_skip:
            # skip unwanted samples
            logger.info(f'{datetime.datetime.now()}: skipping {dir_path} (according to user\'s request)')
            continue
        try:
            with open(os.path.join(dir_path, 'meme.txt')) as f:
                if not first_meme:
                    # skip meme header (6 rows) and add space before next pssm
                    for i in range(6):
                        f.readline()
                merged_meme_f.write(f.read())
                first_meme = False
        except:
            logger.error(f'{datetime.datetime.now()}: no meme file in {dir_path}\n'
                         f'This is its content:\n{os.listdir(dir_path)}')

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')



if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('motif_inference_path', help='A path in which each folder corresponds to a sample and contains a meme file for the motifs in this sample.')
    parser.add_argument('biological_condition', help='A biological condition to merge its samples meme files')
    parser.add_argument('merged_meme_path', help='A path to the output file')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--skip_sample', default='a_weird_str_that_shouldnt_be_a_sample_name_by_any_chance',
                        help='A sample name that should be skipped, e.g., for testing purposes. More than one sample '
                             'name should be separated by commas but no spaces. '
                             'For example: 17b_05,17b_05_test,another_one')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    merge_meme_files(args.motif_inference_path, args.biological_condition, args.merged_meme_path,
                     args.done_file_path, args.skip_sample.split(','), sys.argv)
