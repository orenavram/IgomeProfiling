import os
from sys import argv

path_to_manipulate = argv[1]

for file_name in os.listdir(path_to_manipulate):
    if file_name.endswith('txt'):
        in_path = os.path.join(path_to_manipulate, file_name)
        out_path = os.path.join(path_to_manipulate, file_name[:-3]+'united_peptides.txt')
        records = {}
        with open(in_path) as f:
            for line in f:
                if line.startswith('>'):
                    counts = int(line.rstrip().split()[1])
                else:
                    seq = line.rstrip()
                    records[seq] = records.get(seq, 0) + counts

        sorted_seqs_by_counts = sorted(records, key=records.get, reverse=True)
        result = ''
        for i in range(len(sorted_seqs_by_counts)):
            seq = sorted_seqs_by_counts[i]
            result += f'>P{i+1} {records[seq]}\n{seq}\n'
        with open(out_path, 'w') as f:
            f.write(result)



