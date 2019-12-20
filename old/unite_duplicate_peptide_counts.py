
###############################################################################
# must run only after in pepsurf format, using convert_to_pepsurf.py script!! #
###############################################################################

import os
from sys import argv

path_to_manipulate = argv[1]
path_to_manipulate = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_2/analyses/Avia_TB/random_forest/SICK/Top_11/SEQS'

for file_name in os.listdir(path_to_manipulate):
    if file_name.endswith('txt'):
        in_path = os.path.join(path_to_manipulate, file_name)
        out_path = os.path.join(path_to_manipulate, file_name[:-3]+'united_peptides.txt')
        records = {}
        print(f'Processing {in_path}...')
        with open(in_path) as f:
            for line in f:
                if line.startswith('>'):
                    counts = int(line.rstrip().split()[1])
                else:
                    seq = line.rstrip()
                    if seq in records:
                        print(f'Unifying (at least) second appearance of {seq}')
                    records[seq] = records.get(seq, 0) + counts

        sorted_seqs_by_counts = sorted(records, key=records.get, reverse=True)
        result = ''
        for i in range(len(sorted_seqs_by_counts)):
            seq = sorted_seqs_by_counts[i]
            result += f'>P{i+1} {records[seq]}\n{seq}\n'
        with open(out_path, 'w') as f:
            f.write(result)



