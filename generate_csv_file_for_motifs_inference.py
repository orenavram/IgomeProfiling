import os
from sys import argv
first_phase_output_dir = argv[1]
resuls_path = argv[2]
result=''
for sample_name in os.listdir(first_phase_output_dir):
    sample_path=os.path.join(first_phase_output_dir, sample_name)
    if not os.path.isdir(sample_path):
        continue
    file = ''
    for file_name in os.listdir(sample_path):
        if 'UNIQUE.NORMALIZATION_RpM' in file_name:
            break
    if not file_name:
        raise ValueError
    file_path = os.path.join(sample_path, file_name)
    result += file_path + ',' + sample_name.split('_')[2] + ',' + sample_name + '\n'
with open(resuls_path, 'w') as f:
    f.write(result)