import subprocess
import sys
import os

folder = sys.argv[1]
folder = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_2/analyses/Avia_TB/random_forest/SICK/Top_11/SEQS'
for file in os.listdir(folder):
    in_file_path = os.path.join(folder, file)
    out_file_path = in_file_path[:in_file_path.rindex('.')] + '.txt' # change suffix
    subprocess.call(f'perl /Users/Oren/Dropbox/Projects/gershoni/src/prepare_to_pepsurf.pl {in_file_path} {out_file_path} ALL', shell=True)
