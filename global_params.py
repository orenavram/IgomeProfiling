# modules to load
python = "python/python-anaconda3.6.5"
gcc = "gcc/gcc-8.2.0"
mafft = "mafft/7.123"

# external scripts
pssm_score_peptide_script = "/groups/pupko/orenavr2/gershoni/src/PSSM_score_Peptide_Jan2018.NonVerbose/PSSM_score_Peptide.verbose"
qsub_script = "/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py"

# user input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_path', type=lambda x: x.rstrip('/'), help='A path to the data folder')
parser.add_argument('output_path', type=lambda x: x.rstrip('/'), help='A path to which the results will be written')
parser.add_argument('fastq', help='A path to the fastq file')
parser.add_argument('reads_type', choices=['short', 'long'],
                    help='Is it an old chip (e.g., exp12) or new one (e.g., DP)?')
parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to',
                    choices=['pupkoweb', 'pupkolab', 'pupkotmp',
                             'pupkowebr', 'pupkolabr', 'pupkotmpr',
                             'itaym', 'lilach', 'bioseq', 'bental',
                             'oren.q', 'bioseq20.q'], default='pupkolab')

args = parser.parse_args()

folder_input = args.input_path
folder_output = args.output_path
fastq_file = args.fastq
queue_name = args.q


src_dir = "/groups/pupko/orenavr2/gershoni/src/IgOmeProfiling"
biggest_cluster = 100
biggest_cluster_sec = 400
max_peptide = 1000
gap_threshold = 0.1  # alignment columns with more than gap_threshold proportion of gaps are discarded
split_num = 10
gz = "no"
MistakeAllowed = 1 #default # default == 1
