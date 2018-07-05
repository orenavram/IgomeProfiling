import os
from sys import argv
first_phase_output_dir = argv[1]
resuls_path = argv[2]
result=''
for sample_name in os.listdir(first_phase_output_dir):
    sample_path=os.path.join(first_phase_output_dir, sample_name)
    print('sample_name is: ' + str(sample_name))
    if not os.path.isdir(sample_path):
        print('Skipping sample_name ' + str(sample_name) + ' (not a directory).')
        continue
    file = ''
    file_name_is_legal = False
    for file_name in os.listdir(sample_path):
        if 'UNIQUE.NORMALIZATION_RpM' in file_name:
            file_name_is_legal = True
            break
    if not file_name_is_legal:
        print('sample_name ' + str(sample_name) + ' does not contain UNIQUE.NORMALIZATION_RpM file.')
        raise ValueError
    file_path = os.path.join(sample_path, file_name)
    print('file_path is: ' + str(file_path))
    #if len(sample_name)>8 and sample_name[8] == '_': #e.g., TTACAGCG_Israel_HCV_RNAplus_02_a_Random
    result += file_path + ',' + sample_name.split('_')[2] + ',' + sample_name + '\n'
with open(resuls_path, 'w') as f:
    f.write(result)