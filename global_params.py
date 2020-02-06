import os

is_run_on_cluster = True
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
    is_run_on_cluster = False

local_command_prefix = "bash"

# modules to load
python = "python/python-anaconda3.6.5"
gcc = "gcc/gcc-8.2.0"
mafft = "mafft/7.123"

# external scripts
pssm_score_peptide_script = "/groups/pupko/orenavr2/igomeProfilingPipeline/src/PSSM_score_Peptide/PSSM_score_Peptide"
qsub_script = "/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py"

biggest_cluster = 100
biggest_cluster_sec = 400
max_peptide = 1000
maximal_gap_frequency_allowed_per_column = 0.1  # alignment columns with more than gap_threshold proportion of gaps are discarded
split_num = 10
MistakeAllowed = 1 # default

if __name__ == '__main__':
    # user input
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', type=lambda x: x.rstrip('/'), help='A path to the data folder')
    parser.add_argument('output_path', type=lambda x: x.rstrip('/'), help='A path to which the results will be written')
    parser.add_argument('fastq', help='A path to the fastq file')
    parser.add_argument('barcode2samplename', help='A path to a table with barcode to sample name mapping')
    parser.add_argument('samplename2biologicalcondition',
                        help='A path to a table with sample name to biological condition mapping')

    parser.add_argument('-q', '--queue_name', help='The cluster to which the job(s) will be submitted to',
                        choices=['pupkoweb', 'pupkolab', 'pupkotmp',
                                 'pupkowebr', 'pupkolabr', 'pupkotmpr',
                                 'itaym', 'lilach', 'bioseq', 'bental',
                                 'oren.q', 'bioseq20.q'], default='pupkolab')

    args = parser.parse_args()

    folder_input = args.input_path
    folder_output = args.output_path
    fastq_file = args.fastq
    barcode2samplename_path = args.barcode2samplename
    samplename2biologicalcondition_path = args.samplename2biologicalcondition
    queue_name = args.q

