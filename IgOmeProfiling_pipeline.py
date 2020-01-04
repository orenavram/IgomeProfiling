import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import *

def run_pipeline(fastq_path, barcode2samplename_path, samplename2biologicalcondition_path, analysis_dir, logs_dir,
                 left_construct, right_construct, max_mismatches_allowed, min_sequencing_quality, gz,
                 max_msas_per_sample, max_msas_per_bc,
                 number_of_random_pssms,
                 run_summary_path, error_path, queue, verbose, argv):

    os.makedirs(os.path.split(run_summary_path)[0], exist_ok=True)

    f_run_summary_path = open(run_summary_path, 'w')
    f_run_summary_path.write(' '.join(argv) + '\n\n')
    f_run_summary_path.flush()

    start_time = datetime.datetime.now()

    exp_name = analysis_dir.rstrip('/').split('/')[-2]

    # output folders of the different modules
    first_phase_output_path = os.path.join(analysis_dir, 'reads_filtration')
    second_phase_output_path = os.path.join(analysis_dir, 'motif_inference')
    third_phase_output_path = os.path.join(analysis_dir, 'model_fitting')

    first_phase_done_path = f'{logs_dir}/reads_filtration_done.txt'
    if not os.path.exists(first_phase_done_path):
        os.makedirs(first_phase_output_path, exist_ok=True)
        first_phase_logs_path = os.path.join(logs_dir, 'reads_filtration')
        os.makedirs(first_phase_logs_path, exist_ok=True)

        module_parameters = [fastq_path, first_phase_output_path, first_phase_logs_path,
                             barcode2samplename_path, left_construct, right_construct,
                             max_mismatches_allowed, min_sequencing_quality, first_phase_done_path,
                             '--gz' if gz else '', f'--error_path {error_path}', '-v' if verbose else '']
        cmd = submit_pipeline_step(f'{src_dir}/reads_filtration/module_wraper.py',
                             [module_parameters],
                             logs_dir, f'{exp_name}_reads_filtration',
                             queue, verbose)

        wait_for_results('reads_filtration', logs_dir, num_of_expected_results=1, example_cmd=cmd,
                         error_file_path=error_path, suffix='reads_filtration_done.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping reads filtration. Done file exists at:\n{first_phase_done_path}')

    second_phase_done_path = f'{logs_dir}/motif_inference_done.txt'
    if not os.path.exists(second_phase_done_path):
        os.makedirs(second_phase_output_path, exist_ok=True)
        second_phase_logs_path = os.path.join(logs_dir, 'motif_inference')
        os.makedirs(second_phase_logs_path, exist_ok=True)

        module_parameters = [first_phase_output_path, second_phase_output_path, second_phase_logs_path,
                             samplename2biologicalcondition_path, max_msas_per_sample, max_msas_per_bc,
                             second_phase_done_path,
                             f'--error_path {error_path}', '-v' if verbose else '', f'-q {queue}']
        cmd = submit_pipeline_step(f'{src_dir}/motif_inference/module_wraper.py',
                             [module_parameters],
                             logs_dir, f'{exp_name}_motif_inference',
                             queue, verbose)

        wait_for_results('motif_inference', logs_dir, num_of_expected_results=1, example_cmd=cmd,
                         error_file_path=error_path, suffix='motif_inference_done.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping motif inference. Done file exists at:\n{second_phase_done_path}')

    third_phase_done_path = f'{logs_dir}/model_fitting_done.txt'
    if not os.path.exists(third_phase_done_path):
        os.makedirs(third_phase_output_path, exist_ok=True)
        third_phase_logs_path = os.path.join(logs_dir, 'model_fitting')
        os.makedirs(third_phase_logs_path, exist_ok=True)

        module_parameters = [first_phase_output_path, second_phase_output_path, third_phase_output_path,
                             third_phase_logs_path, samplename2biologicalcondition_path, number_of_random_pssms,
                             third_phase_done_path, f'--error_path {error_path}', '-v' if verbose else '',
                             f'-q {queue}']
        cmd = submit_pipeline_step(f'{src_dir}/model_fitting/module_wraper.py',
                             [module_parameters],
                             logs_dir, f'{exp_name}_model_fitting',
                             queue, verbose)

        wait_for_results('model_fitting', logs_dir, num_of_expected_results=1, example_cmd=cmd,
                         error_file_path=error_path, suffix='model_fitting_done.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping model fitting. Done file exists {third_phase_done_path}')

    end_time = datetime.datetime.now()
    f_run_summary_path.write(f'Total running time: {str(end_time-start_time)[:-3]}')
    f_run_summary_path.close()

    logger.info(f'Started running at {start_time}')
    logger.info(f'Done running at {end_time}')
    logger.info(f'Total running time: {str(end_time-start_time)[:-3]}')
    logger.info('Bye!')



if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}', flush=True)

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('fastq_path', type=str, help='A fastq file to parse')
    parser.add_argument('barcode2samplename_path', type=str, help='A path to the barcode to sample name file')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='A path to the sample name to biological condition file')
    parser.add_argument('analysis_dir', type=str, help='analysis folder')
    parser.add_argument('logs_dir', type=str, help='logs folder')

    # optional parameters for the filtration step
    parser.add_argument('--left_construct', default='CAACGTGGC', help='The (constant) sequence from the LEFT of the random sequence') # in exp12: "CAACGTGGC"
    parser.add_argument('--right_construct', default='GCCT', help='The (constant) sequence from the RIGHT of the random sequence') # in exp12: "GCCT"
    parser.add_argument('--max_mismatches_allowed', type=int, default=1,
                        help='number of mismatches allowed together in both constant sequences')
    parser.add_argument('--min_sequencing_quality', type=int, default=38,
                        help='Minimum average sequencing threshold allowed after filtration'
                             'for more details, see: https://en.wikipedia.org/wiki/Phred_quality_score')
    parser.add_argument('--gz', action='store_true', help='gzip fastq, filtration_log, fna, and faa files')

    # optional parameters for the motif inference
    parser.add_argument('--max_msas_per_sample', default=100, type=int,
                        help='For each sample, align only the biggest $max_msas_per_sample')
    parser.add_argument('--max_msas_per_bc', default=400, type=int,
                        help='For each biological condition, align only the biggest $max_msas_per_bc')

    # optional parameters for the modelling step
    parser.add_argument('--number_of_random_pssms', default=100, type=int, help='Number of pssm permutations')

    # general optional parameters
    parser.add_argument('--run_summary_path', type=str,
                        help='a file in which the running configuration and timing will be written to')
    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('-q', '--queue', default='pupkoweb', type=str, help='a queue to which the jobs will be submitted')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    run_summary_path = args.error_path if args.error_path else os.path.join(args.analysis_dir, 'run_summary_path.txt')
    error_path = args.error_path if args.error_path else os.path.join(args.logs_dir, 'error.txt')

    run_pipeline(args.fastq_path, args.barcode2samplename_path, args.samplename2biologicalcondition_path,
                 args.analysis_dir.rstrip('/'), args.logs_dir.rstrip('/'),
                 args.left_construct, args.right_construct, args.max_mismatches_allowed, args.min_sequencing_quality, True if args.gz else False,
                 args.max_msas_per_sample, args.max_msas_per_bc,
                 args.number_of_random_pssms,
                 run_summary_path, error_path, args.queue, True if args.verbose else False, sys.argv)


'''

pssm_score_peptide = "/groups/pupko/orenavr2/gershoni/src/PSSM_score_Peptide_Jan2018.NonVerbose/PSSM_score_Peptide.verbose"
python = "python/python-anaconda3.6.5"
gcc = "gcc/gcc-8.2.0"
qsub = f"python {src_dir}/q_submitter_power.py"

input_path = args.input_path.rstrip('/')
output_path = args.output_path.rstrip('/')
first_phase_path = f'{output_path}/first_phase_output'
first_phase_logs = f'{first_phase_path}/logs'
motif_inference_path = f'{output_path}/motif_inference'
motif_inference_logs = f'{motif_inference_path}/logs'
classificaion_path = f'{output_path}/model_fitting'
classificaion_logs = f'{classificaion_path}/logs'

os.makedirs(output_path, exist_ok=True)

os.makedirs(first_phase_path, exist_ok=True)
os.makedirs(f'{first_phase_path}/logs', exist_ok=True)

os.makedirs(motif_inference_path, exist_ok=True)
os.makedirs(f'{motif_inference_path}/logs', exist_ok=True)

os.makedirs(classificaion_path, exist_ok=True)
os.makedirs(classificaion_logs, exist_ok=True)
os.makedirs(f'{classificaion_logs}/done', exist_ok=True)
os.makedirs(f'{classificaion_path}/pvalues', exist_ok=True)
os.makedirs(f'{classificaion_path}/hits', exist_ok=True)

configuration_file_path = f'{output_path}/configuration.txt'
logger.info(f'Writing configuration file to: {configuration_file_path}')
#TODO: write configuration file properly
# with open(configuration_file_path, 'w') as f:
#     with open("/groups/pupko/bialik/py/new/global_params.py", 'r') as f2:
#         for l in f2:
#             f.write(l)
#
#
# logger.info(f'{datetime.datetime.now()}: Loading samplename2biologicalcondition file from:\n{samplename2biologicalcondition_path}')
# biological_conditions = set()
# with open(samplename2biologicalcondition_path, 'r') as f:
#     for line in f:
#         bc = line.rstrip().split('\t')[1]
#         biological_conditions.add(bc)
#

# Filter reads step




# Motif inference step


# Classification step

logger.info(f'{datetime.datetime.now()}: Filtering the following fastq file:\n{fastq_file}')
subprocess.call(
    f"python {src_dir}/reads_filtration.py " 
    f"{fastq_file} "
    f"{first_phase_path} "
    f"{barcode2samplename_path} ",
    shell=True)


# TODO: wait for the results

# TODO: need to change it to sent the coutuniqseq_fasta.py as jobs
logger.info(f'{datetime.datetime.now()}: Generating unique_rpm.faa files for every faa file '
            f'(i.e., normalized and without duplications) with and without flanking Cysteine')
for path, dirs, filenames in os.walk(folder_output + "first_phase_output"):
    for filename in filenames:
        if filename.endswith('faa'):
            in_file = os.path.join(path, filename)
            prefix, suffix = os.path.splitext(in_file)
            out_file = f'{prefix}_unique_rpm{suffix}'
            subprocess.call(f"python {src_dir}/count_and_collapse_duplicates.py "
                            f"{in_file} "
                            f"{out_file} "
                            f"--rpm {first_phase_path}/rpm_factors.txt ", shell=True)

            subprocess.call(f"python {src_dir}/remove_cysteine_loop.py "
                            f"{out_file} "  # This is OK (previous out_file is now in_file)
                            f"{prefix}_unique_rpm_cysteine_trimmed{suffix}", shell=True)


#TODO: create a done path
# TODO: wait for the results

logger.info(f'{datetime.datetime.now()}: Inferring motifs...')
subprocess.call(
    f"python /groups/pupko/bialik/py/new/motif_inference.py {folder_input}/metadata/ {folder_output}/first_phase_output/ {folder_output}/motif_inference"
    f" {queue_name} {max_peptide} {biggest_cluster} {biggest_cluster_sec} {maximal_gap_frequency_allowed_per_column}", shell=True)
# print(datetime.datetime.now())


# sequence logo
with open(f'{motif_inference_logs}/weblogo.cmds','w') as f:
    for bc in biological_conditions:
        os.makedirs(f'{motif_inference_path}/{bc}/weblogo', exist_ok=True)
        for alnfile in os.listdir(f'{motif_inference_path}/{bc}/fixed_alignments'):
            aln_file_prefix = os.path.splitext(alnfile)[0]
            f.write(f"module load {python} {gcc}!@#")
            f.write(f'{src_dir}/generate_weblogo.py {motif_inference_path}/{bc}/fixed_alignments/{alnfile} '
                    f'{motif_inference_path}/{bc}/weblogo/{aln_file_prefix}.png!@#\t{aln_file_prefix}\n')
cmd = f'{qsub} {motif_inference_logs}/weblogo.cmds {motif_inference_logs} -q {queue_name}'
subprocess.call(cmd, shell=True)


# print(datetime.datetime.now())
append_write = 'w'
jobs = 0
samples_tup = []
for path, dirname, filename in os.walk(output_path + "/first_phase_output"):
    for f in filename:
        if f.endswith("RpM.fs"):  # 'RpM.TrimmedCysLoop.fas'):
            file = path + '/' + f
            samples_tup.append((file, os.stat(file).st_size))
samples = []
for i in sorted(samples_tup, key=lambda x: x[1], reverse=True):
    samples.append(i[0])

for bc in biological_conditions:
            if "test" in bc:
                continue
            if not os.path.exists(output_path + "/model_fitting/logs/" + bc):
                os.mkdir(output_path + "/model_fitting/logs/" + bc)
            subprocess.call(
                f"python {src_dir}/split_meme_like_motifs_files_and_cutoffs.py {output_path}/motif_inference/{bc} {split_num}", shell=True)
            for i in range(len(os.listdir(output_path + "motif_inference/" + bc + "/pssm"))):
                pssm = f"{output_path}motif_inference/{bc}/pssm/{i}.txt"
                cutoffs = f"{output_path}motif_inference/{bc}/cutoffs/{i}.txt"
                for sample in samples:                    
                    sample_name = sample.split('/')[len(sample.split('/'))-2].replace('_','-')
                    with open(output_path + "/model_fitting/Calc_Pval_for_motifs_vs_samples.S100.cmds", append_write) as f:
                        f.write(f"module load {python} {gcc}!@#{pssm_score_peptide} -pssm {pssm}"
                                f" -pssm_cutoffs {cutoffs} -seq {sample}"
                                f" -out {output_path}model_fitting/logs/{bc}/{bc}"
                                f"_motifs_vs_sample_{sample_name}_{i}.S100.txt -NrandPSSM 100 -CalcPSSM_Pval!@#"
                                f"touch {output_path}/model_fitting/logs/done/{jobs}.done!@#\t{bc}_vs_{sample_name}_S100\n")
                    append_write = 'a'
                    jobs += 1

subprocess.call(f"{qsub} {output_path}/model_fitting/Calc_Pval_for_motifs_vs_samples.S100.cmds {output_path}/model_fitting/logs/ -q {queue_name}", shell=True)
for bc in biological_conditions:
    while True:
        if len(os.listdir(output_path + "/model_fitting/logs/done")) == jobs:
            break
        else:
            print(f"{len(os.listdir(output_path + '/model_fitting/logs/done'))} jobs done out of {jobs}")
            time.sleep(60)
# print(datetime.datetime.now())

print("PSSM_score_Peptide has been finished")

# print(datetime.datetime.now())
for bc in biological_conditions:
    subprocess.call(f"python {src_dir}/pvalues_new.py {output_path}/model_fitting/logs/{bc} "
                    f"{output_path}/model_fitting/ "
                    f"{samplename2biologicalcondition_path}", shell=True)


# print(datetime.datetime.now())

csv_files = []
print("pvalues has been finished")
for path, dirname, files in os.walk(output_path + "model_fitting"):
    for file in files:
        if file.endswith(".csv"):
            csv_files.append(path + '/' + file)

logger.info(f'Applying random forest...')
# print(datetime.datetime.now())
for csv in csv_files:
    subprocess.call(f"python {src_dir}/random_forest.py {csv}", shell=True)
# print(datetime.datetime.now())

# biological_condition = ""
# with open(input_path + "/metadata/samplename2biologicalcondition.txt", 'r') as f:
#     for line in f:
#         line = line.rstrip().split("\t")[1]
#         if line not in biological_condition and "test" not in line:
#             biological_condition += line + ','
# biological_condition = biological_condition.rstrip(',')

for i in ['hits','pvalues']:
    subprocess.call(f"python {src_dir}/merge_biologicalconditions_scores.py {output_path}model_fitting/{i} {biological_condition}", shell=True)



motif_folder = f"{folder_output}/motif_inference/"
with open(folder_output + "motif_inference/logs/weblogo.cmds",'w') as f:
    for folder in os.listdir(motif_folder):
        if folder == "logs":
            continue
        if not os.path.exists(f"{motif_folder}/{folder}/weblogo"):
            os.mkdir(f"{motif_folder}/{folder}/weblogo")
        for file in os.listdir(f"{motif_folder}/{folder}/fixed_alignments/"):
            f.write(f"module load {python}!@#")
            f.write(f"{python_scripts}generate_weblogo.py "
                    f"{motif_folder}/{folder}/fixed_alignments/{file} "
                    f"{motif_folder}/{folder}/weblogo/{file.rstrip('aln')}pdf"
                    f"\t{file.rstrip('.aln')}\n")

cmd = f"{qsub} {motif_inference_path}/logs/weblogo.cmds {motif_inference_path}/logs -q {queue_name}"
#subprocess.call(cmd, shell=True)
'''