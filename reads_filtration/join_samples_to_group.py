import sys


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}', flush=True)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('samples_path_reads', type=str, help='A path to a folder that contain the reads samples')
    parser.add_argument('output_path_group', type=str, help='A path to folder to output the files of the new group')
    parser.add_argument('sample_names', type=str, help='sample names that in the group, the samples seperate by comma, i.e. 17b_01,17_02')
    parser.add_argument('done_file_path', type=str, help='A path to a file that signals that the script finished running successfully')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')
