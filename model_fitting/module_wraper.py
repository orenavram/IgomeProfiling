import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import *
from auxiliaries.stop_machine_aws import stop_machines
from auxiliaries.validation_files import is_input_files_valid 

def repeat_items(list):
    output = []
    for x in list:
        output.append(x)
        output.append(x)
    return output
                      
def build_classifier(first_phase_output_path, motif_inference_output_path,
                     classification_output_path, logs_dir, samplename2biologicalcondition_path, number_of_random_pssms,
                     fitting_done_path, check_files_valid, stop_before_random_forest, run_local_randon_forst_for_bc, num_of_random_configurations_to_sample,
                     number_parallel_random_forest, min_value_error_random_forest,
                     rank_method, tfidf_method, tfidf_factor, shuffles, shuffles_percent, shuffles_digits,
                     cv_num_of_splits, random_forest_seed, random_forest_seed_configurations,
                     stop_machines_flag, type_machines_to_stop, name_machines_to_stop,
                     queue_name, verbose, error_path, use_mapitope, argv):                  

    if check_files_valid and not is_input_files_valid(samplename2biologicalcondition_path=samplename2biologicalcondition_path, barcode2samplename_path='', logger=logger):
        return
                     

    is_pval = rank_method == 'pval'
    os.makedirs(classification_output_path, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    if os.path.exists(fitting_done_path):
        logger.info(f'{datetime.datetime.now()}: skipping model_fitting step ({fitting_done_path} already exists)')
        return

    samplename2biologicalcondition = load_table_to_dict(samplename2biologicalcondition_path, 'Barcode {} belongs to more than one sample_name!!')
    sample_names = sorted(samplename2biologicalcondition)
    biological_conditions = sorted(set(samplename2biologicalcondition.values()))

    for bc in biological_conditions:
        bc_dir_path = os.path.join(classification_output_path, bc)
        os.makedirs(bc_dir_path, exist_ok=True)
        scanning_dir_path = os.path.join(bc_dir_path, 'scanning')
        os.makedirs(scanning_dir_path, exist_ok=True)

    # compute scanning scores (hits and values)
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: scanning peptides vs motifs (hits and values)')
    script_name = 'scan_peptides_vs_motifs.py'
    num_of_expected_results = 0
    num_of_cmds_per_job = 4
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        # get current biological condition (splitted) motifs folder
        memes_path = os.path.join(motif_inference_output_path, bc, 'memes')
        cutoffs_path = os.path.join(motif_inference_output_path, bc, 'cutoffs')

        for file_name in sorted(os.listdir(memes_path)):
            # extract each split of the motifs to scan peptides against it
            meme_file_path = os.path.join(memes_path, file_name)
            cutoffs_file_path = os.path.join(cutoffs_path, file_name)
            for sample_name in sample_names:
                sample_first_phase_output_path = os.path.join(first_phase_output_path, sample_name)
                faa_file_path = get_faa_file_name_from_path(sample_first_phase_output_path, use_mapitope)
                output_path = os.path.join(classification_output_path, bc, 'scanning',
                                           f'{sample_name}_peptides_vs_{bc}_motifs_{os.path.splitext(file_name)[0]}.txt')
                done_path = os.path.join(logs_dir, f'{sample_name}_peptides_vs_{bc}_motifs_{os.path.splitext(file_name)[0]}_done_scan.txt')
                if not os.path.exists(done_path):
                    cmd = [meme_file_path, cutoffs_file_path, faa_file_path, rank_method, str(number_of_random_pssms)]
                    if rank_method == 'shuffles':
                        cmd += ['--shuffles', shuffles]
                    cmd += ['--shuffles_percent', shuffles_percent, '--shuffles_digits', shuffles_digits, output_path, done_path]
                    all_cmds_params.append(cmd)
                else:
                    logger.debug(f'skipping scan as {done_path} found')
                    num_of_expected_results += 1

    if len(all_cmds_params) > 0:
        for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
            current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
            done_path_index = -1 if is_pval else -2
            done_file_name = os.path.split(current_batch[0][done_path_index])[-1]
            name_tokens = done_file_name.split('_peptides_vs_')
            logger.info(name_tokens)
            sample_name = name_tokens[0]
            bc = name_tokens[1].split('_motifs_')[0]
            split_num = name_tokens[1].split('_motifs_')[1].split('_done_')[0]
            assert sample_name in sample_names, f'Sample {sample_name} not in sample names list:\n{sample_names}'
            assert bc in biological_conditions, f'Biological condition {bc} not in bc names list:\n{biological_conditions}'
            cmd = submit_pipeline_step(f'{src_dir}/model_fitting/{script_name}', current_batch,
                                    logs_dir, f'{sample_name}_vs_{bc}_scan_{split_num}', queue_name, verbose)
            num_of_expected_results += len(current_batch)

        wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path=error_path, suffix='_done_scan.txt')
    else:
        logger.info(f'Skipping scanning peptides vs motifs (hits and values), all scans found')


    # aggregate scanning scores (hits and values)
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: aggregating scores')
    script_name = 'merge_pvalues.py'
    num_of_expected_results = 0
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        meme_path = os.path.join(motif_inference_output_path, bc, 'meme.txt')
        scanning_dir_path = os.path.join(classification_output_path, bc, 'scanning')
        done_path = os.path.join(logs_dir, f'{bc}_done_aggregate_scores.txt')
        if not os.path.exists(done_path):
            if not is_pval:
                cmds = ['--memes', meme_path,
                        '--bc', bc,
                        '--sam2bc', samplename2biologicalcondition_path,
                        '--scan', scanning_dir_path,
                        '--output', classification_output_path,
                        '--method', tfidf_method,
                        '--factor', str(tfidf_factor),
                        '--done', done_path]
                if rank_method == 'shuffles':
                    cmds.append('--rank')
                all_cmds_params.append(cmds)
            else:
                aggregated_values_path = os.path.join(classification_output_path, bc, f'{bc}_values.csv')
                aggregated_hits_path = os.path.join(classification_output_path, bc, f'{bc}_hits.csv')
                all_cmds_params.append([meme_path, scanning_dir_path, bc, aggregated_values_path,
                                        aggregated_hits_path, samplename2biologicalcondition_path, done_path])
        else:
            logger.debug(f'skipping score aggregation as {done_path} found')
            num_of_expected_results += 1

    if len(all_cmds_params) > 0:
        executable = 'python' if is_pval else None
        script_path = f'{src_dir}/model_fitting/{script_name}' if is_pval else f'{src_dir}/tfidf/tfidf'
        for cmds_params, bc in zip(all_cmds_params, biological_conditions):
            cmd = submit_pipeline_step(script_path,[cmds_params],
                                logs_dir, f'{bc}_aggregate_scores',
                                queue_name, verbose, executable=executable)
            num_of_expected_results += 1  # a single job for each biological condition

        wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path=error_path, suffix='_done_aggregate_scores.txt')
    else:
        logger.info(f'skipping aggregating scores, all scores found')

    
    # fitting a random forest model (hits and values)
    if not stop_before_random_forest:
        logger.info('_'*100)
        logger.info(f'{datetime.datetime.now()}: fitting model')
        script_name = 'random_forest.py'
        num_of_expected_results = 0
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for bc in biological_conditions:
            aggregated_values_path = os.path.join(classification_output_path, bc, f'{bc}_values.csv')
            pvalues_done_path = os.path.join(logs_dir, f'{bc}_values_done_fitting.txt')
            aggregated_hits_path = os.path.join(classification_output_path, bc, f'{bc}_hits.csv')
            hits_done_path = os.path.join(logs_dir, f'{bc}_hits_done_fitting.txt')
            
            value_cmd = [aggregated_values_path, pvalues_done_path, logs_dir, error_path, f'--num_of_configurations_to_sample {num_of_random_configurations_to_sample}', f'--cv_num_of_splits {cv_num_of_splits}',
                    f'--number_parallel_random_forest {number_parallel_random_forest}', f'--min_value_error_random_forest {min_value_error_random_forest}', f'--seed {random_forest_seed}',
                    f'--random_forest_seed {random_forest_seed_configurations}', f'--rank_method {rank_method}', f'--queue {queue_name}']
            hits_cmd = [aggregated_hits_path, hits_done_path, logs_dir, error_path, f'--num_of_configurations_to_sample {num_of_random_configurations_to_sample}', f'--cv_num_of_splits {cv_num_of_splits}',
                    f'--number_parallel_random_forest {number_parallel_random_forest}', f'--min_value_error_random_forest {min_value_error_random_forest}', f'--seed {random_forest_seed}',
                    f'--random_forest_seed {random_forest_seed_configurations}', '--rank_method hits', f'--queue {queue_name}']
            if not os.path.exists(pvalues_done_path):
                all_cmds_params.append(value_cmd)
            else:
                logger.debug(f'Skipping fitting as {pvalues_done_path} found')
                num_of_expected_results += 1
            
            if not os.path.exists(hits_done_path):
                all_cmds_params.append(hits_cmd)
            else:
                logger.debug(f'Skipping fitting as {hits_done_path} found')
                num_of_expected_results += 1

        if len(all_cmds_params) > 0:
            doubled_bc = repeat_items(biological_conditions)
            for cmds_params, bc in zip(all_cmds_params, doubled_bc):
                cmd = ''
                if run_local_randon_forst_for_bc:
                    cmd = run_step_locally(f'{src_dir}/model_fitting/{script_name}',
                                    [cmds_params],
                                    logs_dir, f'{bc}_model',
                                    queue_name, verbose)
                else:
                    cmd = submit_pipeline_step(f'{src_dir}/model_fitting/{script_name}',
                                    [cmds_params],
                                    logs_dir, f'{bc}_model',
                                    queue_name, verbose)
                num_of_expected_results += 1  # a single job for each biological condition

            wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                            error_file_path=error_path, suffix='_done_fitting.txt')
        else:
            logger.info(f'Skipping fitting, all found')
    else:
        logger.info(f'Skipping fitting, stop before random forest')

    if stop_machines_flag:
        stop_machines(type_machines_to_stop, name_machines_to_stop, logger)

    # TODO: fix this bug with a GENERAL WRAPPER done_path
    # wait_for_results(script_name, num_of_expected_results)
    with open(fitting_done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


def get_faa_file_name_from_path(path, use_mapitope):
    for file_name in os.listdir(path):
        if file_name.endswith('faa') and 'unique' not in file_name and ('mapitope' in file_name) == use_mapitope:
            file_name = file_name
            break
    return os.path.join(path, file_name)


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('parsed_fastq_results', type=str, help='A path in which each subfolder corresponds to a samplename and contains a collapsed faa file')
    parser.add_argument('motif_inference_results', type=str, help='A path in which there is a subfolder for each bc '
                                                                  'In each such subfolder there is a memes folder and a cutoffs folder')
    parser.add_argument('classification_output_path', type=str, help='output folder')
    parser.add_argument('logs_dir', type=str, help='logs folder')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='A path to the sample name to biological condition file')
    parser.add_argument('number_of_random_pssms', default=100, type=int, help='Number of pssm permutations')
    parser.add_argument('done_file_path', type=str, help='A path to a file that signals that the module finished running successfully')
    
    parser.add_argument('--check_files_valid', action='store_true', help='Need to check the validation of the files (samplename2biologicalcondition_path / barcode2samplenaem)')
    parser.add_argument('--stop_before_random_forest', action='store_true', help='A boolean flag for mark if we need to run the random forest')
    parser.add_argument('--run_local_randon_forst_for_bc', action='store_true', help='Run the BC random forest local and not parallel')
    parser.add_argument('--num_of_random_configurations_to_sample', default=100, type=int, help='How many random configurations of hyperparameters should be sampled?')
    parser.add_argument('--number_parallel_random_forest', default=20, type=int, help='How many random forest configurations to run in parallel')
    parser.add_argument('--min_value_error_random_forest', default=0.0, type=float, help='A random forest model error value for convergence allowing to stop early')
    parser.add_argument('--rank_method', choices=['pval', 'tfidf', 'shuffles'], default='pval', help='Motifs ranking method')
    parser.add_argument('--tfidf_method', choices=['boolean', 'terms', 'log', 'augmented'], default='boolean', help='TF-IDF method')
    parser.add_argument('--tfidf_factor', type=float, default=0.5, help='TF-IDF augmented method factor (0-1)')
    parser.add_argument('--shuffles', default=5, type=int, help='Number of controlled shuffles permutations')
    parser.add_argument('--shuffles_percent', default=0.2, type=float, help='Percent from shuffle with greatest number of hits (0-1)')
    parser.add_argument('--shuffles_digits', default=2, type=int, help='Number of digits after the point to print in scanning files')
    parser.add_argument('--cv_num_of_splits', default=2, type=int, help='How folds should be in the cross validation process? (use 0 for leave one out)')
    parser.add_argument('--seed_random_forest', default=42, type=int, help='Seed number for reconstructing experiments')
    parser.add_argument('--random_forest_seed_configurations', default=123, type=int, help='Random seed value for generating random forest configurations')
    parser.add_argument('--stop_machines', action='store_true', help='Turn off the machines in AWS at the end of the running')
    parser.add_argument('--type_machines_to_stop', default='', type=str, help='Type of machines to stop, separated by comma. Empty value means all machines. Example: t2.2xlarge,m5a.24xlarge')
    parser.add_argument('--name_machines_to_stop', default='', type=str, help='Names (patterns) of machines to stop, separated by comma. Empty value means all machines. Example: worker*')
    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('-q', '--queue', default='pupkoweb', type=str, help='a queue to which the jobs will be submitted')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    parser.add_argument('-m', '--mapitope', action='store_true', help='use mapitope encoding')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger('main')

    error_path = args.error_path if args.error_path else os.path.join(args.parsed_fastq_results, 'error.txt')

    build_classifier(args.parsed_fastq_results, args.motif_inference_results, args.classification_output_path,
                     args.logs_dir, args.samplename2biologicalcondition_path,  args.number_of_random_pssms, args.done_file_path,
                     args.check_files_valid, args.stop_before_random_forest, args.run_local_randon_forst_for_bc, args.num_of_random_configurations_to_sample, 
                     args.number_parallel_random_forest, args.min_value_error_random_forest, args.rank_method,
                     args.tfidf_method, args.tfidf_factor, args.shuffles, args.shuffles_percent, args.shuffles_digits,
                     args.cv_num_of_splits, args.seed_random_forest, args.random_forest_seed_configurations, 
                     args.stop_machines, args.type_machines_to_stop, args.name_machines_to_stop,
                     args.queue, args.verbose, error_path, args.mapitope, sys.argv)
                      