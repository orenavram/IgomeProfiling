from numpy import append
from stat import FILE_ATTRIBUTE_INTEGRITY_STREAM
import sys
import os

def get_new_meme_lines(input_file):
    meme_file = open(input_file, 'r')
    lines = meme_file.readlines()
    new_lines = []
    flag_letter_probability = False
    flag_Background_letter = False
    for line in lines:
        if line.startswith('MOTIF'):
            new_lines.append('\n') 
            flag_Background_letter = False
            new_lines.append(f'MOTIF {line.split()[1]}\n')
            continue
        if flag_letter_probability:
            new_lines.append(line)
            continue
        if flag_Background_letter:
            line_strip = line.strip('\n')
            new_lines.append(line_strip)
            continue
        if line.startswith('Background'):
            flag_Background_letter = True
        elif line.startswith('letter-probability'):
            new_lines.append(line)
            flag_letter_probability = True
        else:
            continue
    return new_lines


def adjustment_meme_file(input_path_meme_files, output_path):
    num_file = 0
    for name_file in os.listdir(input_path_meme_files):
        input_file_path = os.path.join(input_path_meme_files, name_file)
        new_meme_lines = get_new_meme_lines(input_file_path)
        file_name_meme = ''
        if num_file < 9:
            file_name_meme = f'0{num_file}'
        else:
            file_name_meme = f'{num_file}'    
        output_file = open(os.path.join(output_path, f'{file_name_meme}.txt'), 'w')
        output_file.write(f'MEME version 4\n\n'
                          f'ALPHABET= ACDEFGHIKLMNPQRSTVWY\n\n'
                          f'Background letter frequencies\n')
        for line in new_meme_lines:
            output_file.write(line)
        num_file += 1


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path_meme_files', type=str, help='A path to a meme files')
    parser.add_argument('output_path', type=str, help='Path to put the output meme files.')
    args = parser.parse_args()
    
    adjustment_meme_file(args.input_path_meme_files, args.output_path)
    