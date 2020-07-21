# import datetime
import sys
import logging
import os
logger = logging.getLogger('main')


if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

# needs $src_dir in path
from auxiliaries.pipeline_auxiliaries import load_table_to_dict

mapitope = {
    'R': 'B',
    'K': 'B',
    'E': 'J',
    'D': 'J',
    'S': 'O',
    'T': 'O',
    'I': 'U',
    'L': 'U',
    'V': 'U',
    'Q': 'X',
    'N': 'X',
    'W': 'Z',
    'F': 'Z',
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'H': 'H',
    'M': 'M',
    'P': 'P',
    'Y': 'Y'
}

def convert_to_mapitope(sequence):
    return ''.join([mapitope[c] for c in sequence])


def convert_faa_file_by_maptope(input_file_path, output_file_path):
    with open(input_file_path,'r') as input_file, open(output_file_path,'w+') as output_file:
        input_lines = input_file.readlines()
        for line in input_lines:
            if line[0] == '>':
                output_file.write(line)  
            else:
                output_file.write(convert_to_mapitope(line[:-1])+'\n')

def mapitope_convert_all_samples(reads_filtration_path, barcode2samplename_path, done_file_path, argv='no_argv'):
    barcode2samplename = load_table_to_dict(barcode2samplename_path, 'Barcode {} belongs to more than one sample!!')
    sampleList = set(barcode2samplename.values())
    for sample_name in sampleList:
        input_path = f'{reads_filtration_path}/{sample_name}/{sample_name}.faa'
        output_path = f'{reads_filtration_path}/{sample_name}/{sample_name}_mapitope.faa'
        convert_faa_file_by_maptope(input_path,output_path) 

    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('reads_filtration_path', help='path to reads_filtration folder')
    parser.add_argument('barcode2samplename_path', help='path to barcode2sample file')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    mapitope_convert_all_samples(args.reads_filtration_path, args.barcode2samplename_path, args.done_file_path,sys.argv)