import logging
import os
import sys


def fix_labels(input_path, text_to_catch_in_file, text_to_catch_in_path):

    files_to_manipulate = []
    for path, dirs, files in os.walk(input_path):
        for file in files:
            if text_to_catch_in_file not in file:
                continue
            if text_to_catch_in_path and text_to_catch_in_path not in path.split('/'):
                continue
            files_to_manipulate.append([os.path.join(path, file),
                                        os.path.split(path)[-1],
                                        os.path.join(path, file.replace('pvalues', 'pvalues_fixed'))])

    print('\n'.join(f'{x}\t{y}\n{z}' for x,y,z in files_to_manipulate))
    print('*'*70)

    for file, bc, out in files_to_manipulate:
        with open(file) as f:
            result = f.readline()  # header
            for line in f:
                tokens = line.split(',')
                tokens[1] = bc if bc in tokens[0] else 'other'
                result += ','.join(tokens)
        print(f'Writing results to:\n{out}')
        with open(out, 'w') as f:
            f.write(result)


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', type=lambda x: x.rstrip('/'),
                        help='a path in which pvalues files can be found (somewhere)')
    parser.add_argument('biological_conditions_as_str', help='bcs that should be handled (separated by commas)')
    parser.add_argument('text_to_catch_in_file', default='',
                        help='filters out all files without this text in their name')
    parser.add_argument('--text_to_catch_in_path', default='',
                        help='filters out all paths without *exactly* this text as one of the folder along the path')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    # fix_labels(args.input_path, args.biological_conditions_as_str.split(','), args.text_to_catch_in_file, args.text_to_catch_in_path)
    fix_labels(args.input_path, args.text_to_catch_in_file, args.text_to_catch_in_path)


