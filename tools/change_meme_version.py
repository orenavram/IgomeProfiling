import sys
import os

from sklearn.decomposition import FactorAnalysis
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty


def write_to_file(meme_file_path, lines):
    output_file = open(meme_file_path, 'w')
    output_file.write(f'MEME version 4\n\n')
    for line in lines:
        output_file.write(line)


def split_to_files(meme_file_path, splitted_meme_dir, motifs_per_file):
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


def get_num_lines_pssm(line):
    return int(line.split()[5])


def adjustment_meme_file(input_path_meme_file, output_path, motifs_per_file):
    lines = []
    flag_background = False
    meme_file = open(input_path_meme_file, 'r')
    meme_file_lines = meme_file.readlines()
    count_lines = 0
    for num_line,line in enumerate(meme_file_lines):
        if num_line < count_lines and flag_background:
            lines.append(line.strip())
        elif 'ALPHABET' in line:
            lines.append(f'{line}\n')
        elif 'Background letter frequencies' in line:
            lines.append(f'Background letter frequencies\n')
            flag_background = True
            count_lines = num_line + 3
        elif 'position-specific probability matrix' in line:
            flag_background = False
            lines.append(f'\nMOTIF {line.split()[1]}\n')
            line = meme_file_lines[num_line + 2 ]
            num_lines_to_read  = get_num_lines_pssm(line)
            lines.append(line)
            pssm_lines = meme_file_lines[num_line + 3 : num_line + 3+ num_lines_to_read]
            for line in pssm_lines:
                lines.append(f'{line.lstrip()}')
            lines.append('\n')

    meme_file_path = os.path.join(output_path, 'meme.txt')
    write_to_file(meme_file_path, lines)
    if meme_file_path:
        verify_file_is_not_empty(meme_file_path)
    splitted_meme_dir = os.path.join(os.path.split(meme_file_path)[0], 'memes')
    os.makedirs(splitted_meme_dir, exist_ok=True)
    split_to_files(meme_file_path, splitted_meme_dir, motifs_per_file)


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path_meme_file', type=str, help='A path to a meme file that need to change the version')
    parser.add_argument('output_path', type=str, help='Path to put the output meme files.')
    parser.add_argument('--motifs_per_file', type=int, default=5, help='How many pssm will be in one meme file')
    args = parser.parse_args()
    
    adjustment_meme_file(args.input_path_meme_file, args.output_path, args.motifs_per_file)
    