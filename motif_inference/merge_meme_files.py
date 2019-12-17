import datetime
import os
import logging
logger = logging.getLogger('main')


def merge_meme_files(motif_inference_path, biological_condition, meme_file_name, samples_to_skip):
    """
    :param motif_inference_path: A path in which each folder corresponds to a sample and contains a meme file for the
    motifs in this sample.
    :param biological_condition: A biological condition to merge its samples meme files
    :param meme_file_name: An optional file name for the output file
    :return: A merged meme file at $motif_inference_path/$biological_condition/$meme_file_name
    """
    logger.info(f'{datetime.datetime.now()}: merging meme files of {biological_condition}')

    merged_meme_f = open(os.path.join(motif_inference_path, biological_condition, meme_file_name), 'w')
    first_meme = True
    for dir_name in sorted(os.listdir(motif_inference_path)):
        dir_path = os.path.join(motif_inference_path, dir_name)
        if not os.path.isdir(dir_path) or biological_condition not in dir_name:
            # skip file or folders of non-related biological condition
            continue
        if dir_name in samples_to_skip:
            # skip unwanted samples
            logger.info(f'{datetime.datetime.now()}: skipping analysis for {dir_path}')
            continue
        try:
            with open(os.path.join(dir_path, 'meme.txt')) as f:
                if not first_meme:
                    # skip meme header (6 rows) and add space before next pssm
                    for i in range(6):
                        f.readline()
                    # merged_meme_f.write('\n\n')
                merged_meme_f.write(f.read())
                merged_meme_f.write('\n\n')
                first_meme = False
        except:
            logger.error(f'{datetime.datetime.now()}: no meme file in {dir_path}\n'
                         f'This is its content:\n{os.listdir(dir_path)}')


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('motif_inference_path', help='A path in which each folder corresponds to a sample and contains a meme file for the motifs in this sample.')
    parser.add_argument('biological_condition', help='A biological condition to merge its samples meme files')
    parser.add_argument('--merged_meme_file_name', default='merged_meme.txt',
                        help='A merged meme file at $motif_inference_path/$biological_condition/$meme_file_name')
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

    merge_meme_files(args.motif_inference_path, args.biological_condition,
                     args.merged_meme_file_name, args.skip_sample.split(','))
