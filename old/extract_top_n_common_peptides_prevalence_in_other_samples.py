import subprocess
import sys
import os


def extract_top_n_peptides(path, n):
    result = []
    with open(path) as f:
        i = 1
        for line in f:
            if not line.startswith('>'):
                i += 1
                result.append(line.rstrip())
                if i>n:
                    break
    return result


def extract_peptides_records_in_other_sample(top_n_peptides, file_path):
    result = ''
    with open(file_path) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            peptide = f.readline().strip()
            if peptide in top_n_peptides:
                result += f'{header}\n{peptide}\n'
    return result


in_path = sys.argv[1]
in_folder = sys.argv[2]
out_folder = sys.argv[3]
n = int(sys.argv[4])

top_n_peptides = extract_top_n_peptides(in_path, n)

if not os.path.exists(out_folder):
    os.makedirs(out_folder)

for file in os.listdir(in_folder):
    if file.endswith('txt'):
        file_path = os.path.join(in_folder, file)
        peptides_fasta = extract_peptides_records_in_other_sample(top_n_peptides, file_path)
        with open(os.path.join(out_folder, file), 'w') as f:
            f.write(peptides_fasta)

