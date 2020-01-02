import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import *


def align_clean_pssm_weblogo(folder_names_to_handle, max_clusters_to_align,
                             motif_inference_output_path, logs_dir, error_path, queue_name, verbose):
    # For each sample, align each cluster
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: aligning clusters in each sample')
    script_name = 'align_sequences.py'
    num_of_expected_results = 0
    msas_paths = []  # keep all msas' paths for the next step
    num_of_cmds_per_job = 33
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    # max_number_of_leading_zeros = len(str(max_clusters_to_align))
    logger.info(f'folder_names_to_handle:\n{folder_names_to_handle}')
    for folder in folder_names_to_handle:
        path = os.path.join(motif_inference_output_path, folder, 'unaligned_sequences')
        sample_motifs_dir = os.path.split(path)[0]
        sample_name = os.path.split(sample_motifs_dir)[-1]
        assert sample_name in folder_names_to_handle, f'Sample {sample_name} not in folder names list:\n{folder_names_to_handle}'
        aligned_sequences_path = path.replace('unaligned_sequences', 'aligned_sequences')
        msas_paths.append(aligned_sequences_path)
        os.makedirs(aligned_sequences_path, exist_ok=True)
        for i, faa_filename in enumerate(sorted(os.listdir(path))):  # sorted by clusters rank
            if i == max_clusters_to_align:
                break
            unaligned_cluster_path = os.path.join(path, faa_filename)
            cluster_rank = get_cluster_rank_from(faa_filename)
            # cluster_rank = str(get_cluster_rank_from(faa_filename)).zfill(max_number_of_leading_zeros)
            # tokens = faa_filename.split('_')
            # tokens[tokens.index('clusterRank')+1] = cluster_rank  # update cluster rank with respect to biological condition (and not with respect to sample)
            # aligned_cluster_path = os.path.join(aligned_sequences_path, '_'.join(tokens))
            aligned_cluster_path = os.path.join(aligned_sequences_path, faa_filename)
            done_path = f'{logs_dir}/04_{sample_name}_{cluster_rank}_done_msa.txt'
            all_cmds_params.append([unaligned_cluster_path, aligned_cluster_path, done_path])

    for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
        current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
        sample_name = current_batch[0][1].split('/')[-3]
        assert sample_name in folder_names_to_handle, f'Sample {sample_name} not in folder names list:\n{folder_names_to_handle}'
        cluster_rank = get_cluster_rank_from(current_batch[-1][1])
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             current_batch,
                             logs_dir, f'{sample_name}_{cluster_rank}_msa',
                             queue_name, verbose)

        num_of_expected_results += len(current_batch)


    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='msa.txt')

    # For each sample, clean alignments from gappy columns
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: cleaning alignments from gappy columns')
    script_name = 'remove_gappy_columns.py'
    num_of_expected_results = 0
    cleaned_msas_paths = []  # keep all cleaned msas' paths for the next step
    num_of_cmds_per_job = 50
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    # done_files_list = []
    for msas_path in msas_paths:
        sample_motifs_dir = os.path.split(msas_path)[0]
        sample_name = os.path.split(sample_motifs_dir)[-1]
        assert sample_name in folder_names_to_handle, f'Sample {sample_name} not in folder names list:\n{folder_names_to_handle}'
        cleaned_msas_path = os.path.join(sample_motifs_dir, 'cleaned_aligned_sequences')
        cleaned_msas_paths.append(cleaned_msas_path)
        os.makedirs(cleaned_msas_path, exist_ok=True)
        for msa_name in sorted(os.listdir(msas_path)):  # sorted by clusters rank
            msa_path = os.path.join(msas_path, msa_name)
            cleaned_msa_path = os.path.join(cleaned_msas_path, msa_name)
            done_path = f'{logs_dir}/05_{sample_name}_{msa_name}_done_cleaning_msa.txt'
            # done_files_list.append(done_path)
            all_cmds_params.append([msa_path, cleaned_msa_path, done_path])
    for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
        current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
        sample_name = current_batch[0][1].split('/')[-3]
        assert sample_name in folder_names_to_handle, f'Sample {sample_name} not in folder names list:\n{folder_names_to_handle}'
        cluster_rank = get_cluster_rank_from(current_batch[-1][0])
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             current_batch,
                             logs_dir, f'{sample_name}_{cluster_rank}_clean',
                             queue_name, verbose)
        num_of_expected_results += len(current_batch)

    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='_done_cleaning_msa.txt') #, done_files_list=done_files_list)

    # For each sample, generate a meme file with a corresponding pssm for each alignment
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: generating meme files for each sample from cleaned alignments')
    script_name = 'create_meme.py'
    num_of_expected_memes = 0
    num_of_cmds_per_job = 4
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for cleaned_msas_path in cleaned_msas_paths:
        sample_motifs_dir = os.path.split(cleaned_msas_path)[0]
        sample_name = os.path.split(sample_motifs_dir)[-1]
        assert sample_name in folder_names_to_handle, f'Sample {sample_name} not in folder names list:\n{folder_names_to_handle}'
        meme_path = os.path.join(sample_motifs_dir, 'meme.txt')
        done_path = f'{logs_dir}/06_{sample_name}_done_meme.txt'
        all_cmds_params.append([cleaned_msas_path, meme_path, done_path])
    for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
        current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             current_batch,
                             logs_dir, f'{i//num_of_cmds_per_job}_meme',
                             queue_name, verbose)
        num_of_expected_memes += len(current_batch)

    # instead of waiting here, submit the weblogos first..

    # For each cleaned msa, generate a web logo. No need to wait with the analysis.
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: generating weblogos for each cleaned alignment')
    script_name = 'generate_weblogo.py'
    num_of_cmds_per_job = 100
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for cleaned_msas_path in cleaned_msas_paths:
        weblogos_path = cleaned_msas_path.replace('cleaned_aligned_sequences', 'weblogos')
        os.makedirs(weblogos_path, exist_ok=True)
        for msa_name in os.listdir(cleaned_msas_path):
            msa_path = os.path.join(cleaned_msas_path, msa_name)
            weblogo_prefix_path = os.path.join(weblogos_path, os.path.splitext(msa_name)[0])
            all_cmds_params.append([msa_path, weblogo_prefix_path])

    for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
        current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             current_batch,
                             logs_dir, f'weblogo_{i}th_batch',
                             queue_name, verbose)

    # wait for the memes!! (previous logical block!)
    wait_for_results(script_name, logs_dir, num_of_expected_memes,
                     error_file_path=error_path, suffix='_done_meme.txt')


def infer_motifs(first_phase_output_path, max_msas_per_sample, max_msas_per_bc,
                 motif_inference_output_path, logs_dir, samplename2biologicalcondition_path,
                 motif_inference_done_path, queue_name, verbose, error_path, argv):

    os.makedirs(motif_inference_output_path, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    samplename2biologicalcondition = load_table_to_dict(samplename2biologicalcondition_path, 'Barcode {} belongs to more than one sample!!')
    sample_names = sorted(samplename2biologicalcondition)
    biological_conditions = sorted(set(samplename2biologicalcondition.values()))

    if os.path.exists(motif_inference_done_path):
        logger.info(f'{datetime.datetime.now()}: skipping motif_inference step ({motif_inference_done_path} already exists)')
        return

    # Make sure all sequences are in upper case letters (for example, no need to differentiate between Q and q)
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: upper casing all sequences in the faa files')
    script_name = 'upper_case_sequences.py'
    num_of_expected_results = 0
    upper_faa_paths = [] # keep all faas' paths for the next step
    for sample_name in sample_names:
        dir_path = os.path.join(first_phase_output_path, sample_name)
        assert(os.path.exists(dir_path))
        if not os.path.isdir(dir_path):
            # skip files or folders of non-related biological condition
            continue
        faa_filename = f'{sample_name}_unique_rpm.faa'
        in_faa_path = os.path.join(first_phase_output_path, sample_name, faa_filename)
        out_faa_dir = os.path.join(motif_inference_output_path, sample_name)
        os.makedirs(out_faa_dir, exist_ok=True)
        out_faa_path = os.path.join(out_faa_dir, f'{sample_name}_upper{faa_filename.split(sample_name)[-1]}')
        upper_faa_paths.append(out_faa_path)
        done_path = f'{logs_dir}/01_{sample_name}_done_uppering.txt'
        fetch_cmd(f'{src_dir}/motif_inference/{script_name}',
                  [in_faa_path, out_faa_path, done_path], verbose, error_path)
        num_of_expected_results += 1

    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='uppering.txt')

    # Clustering sequences within each sample
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: clustering sequences in each sample')
    script_name = 'cluster.py'
    num_of_expected_results = 0
    clstr_paths = [] # keep all clusters' paths for the next step
    num_of_cmds_per_job = 1
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for upper_faa_path in upper_faa_paths:
        faa_dir, faa_filename = os.path.split(upper_faa_path)
        sample_name = os.path.split(faa_dir)[-1]
        assert sample_name in sample_names, f'Sample {sample_name} not in sample names list:\n{sample_names}'
        output_prefix = os.path.join(faa_dir, sample_name)
        clstr_paths.append(f'{output_prefix}.clstr')
        done_path = f'{logs_dir}/02_{sample_name}_done_clustering.txt'
        all_cmds_params.append([upper_faa_path, output_prefix, done_path])

    for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
        current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
        sample_name = os.path.split(current_batch[0][1])[-1]
        assert sample_name in sample_names, f'Sample {sample_name} not in sample names list:\n{sample_names}'
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             current_batch,
                             logs_dir, f'{sample_name}_cluster', queue_name, verbose)
        num_of_expected_results += len(current_batch)

    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='clustering.txt')

    # For each sample, split the faa file to the clusters inferred in the previous step
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: extracting clusters sequences for each sample')
    script_name = 'extract_clusters_sequences.py'
    num_of_expected_results = 0
    unaligned_clusters_folders = [] # keep all sequences' paths for the next step
    num_of_cmds_per_job = 1
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for upper_faa_path, clstr_path in zip(upper_faa_paths, clstr_paths):
        faa_dir, faa_filename = os.path.split(upper_faa_path)
        sample_name = os.path.split(faa_dir)[-1]
        assert sample_name in sample_names, f'Sample {sample_name} not in sample names list:\n{sample_names}'
        clusters_sequences_path = os.path.join(faa_dir, 'unaligned_sequences')
        unaligned_clusters_folders.append(clusters_sequences_path)
        os.makedirs(clusters_sequences_path, exist_ok=True)
        done_path = f'{logs_dir}/03_{sample_name}_done_extracting_sequences.txt'
        all_cmds_params.append([upper_faa_path, clstr_path, clusters_sequences_path, done_path,
                                f'--file_prefix {sample_name}'])

    for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
        current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
        clusters_sequences_path = current_batch[0][2]
        sample_name = clusters_sequences_path.split('/')[-2]
        assert sample_name in sample_names, f'Sample {sample_name} not in sample names list:\n{sample_names}'
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             current_batch,
                             logs_dir, f'{sample_name}_extracting_sequences', queue_name, verbose)
        num_of_expected_results += len(current_batch)

    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='extracting_sequences.txt')

    # 3 steps together!! align each cluster; clean each alignment; calculate pssm for each alignment
    align_clean_pssm_weblogo(sample_names, max_msas_per_sample,
                             motif_inference_output_path, logs_dir, error_path, queue_name, verbose)

    # Merge memes of the same biological condition
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: merging meme files for each of the following biological conditions\n'
                f'{biological_conditions}')
    script_name = 'merge_meme_files.py'
    num_of_expected_results = 0
    biological_condition_memes = []
    all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        done_path = f'{logs_dir}/07_{bc}_done_meme_merge.txt'
        bc_folder = os.path.join(motif_inference_output_path, bc)
        os.makedirs(bc_folder, exist_ok=True)
        output_path = os.path.join(bc_folder, 'merged_meme.txt')
        biological_condition_memes.append(output_path)
        all_cmds_params.append([motif_inference_output_path, bc, output_path, done_path])

    for cmds_params, bc in zip(all_cmds_params, biological_conditions):
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             [cmds_params],
                             logs_dir, f'{bc}_merge_meme',
                             queue_name, verbose)
        num_of_expected_results += 1  # a single job for each biological condition

    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='_done_meme_merge.txt')


    # Unite motifs based on their correlation using UnitePSSMs.cpp
    # TODO: verify how exactly this is done
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: detecting pssms to unite for the following biological conditions\n'
                f'{biological_conditions}')
    script_name = 'unite_motifs_of_biological_condition.py'
    num_of_expected_results = 0
    all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for merged_meme_path, bc in zip(biological_condition_memes, biological_conditions):
        relevant_samples = ','.join([sample for sample in samplename2biologicalcondition if samplename2biologicalcondition[sample] == bc])
        output_path = os.path.split(merged_meme_path)[0]
        done_path = f'{logs_dir}/08_{bc}_done_detecting_similar_pssms.txt'
        all_cmds_params.append([motif_inference_output_path, merged_meme_path, bc, relevant_samples, output_path, done_path])

    for cmds_params, bc in zip(all_cmds_params, biological_conditions):
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             [cmds_params],
                             logs_dir, f'{bc}_detect_similar_pssms',
                             queue_name, verbose)
        num_of_expected_results += 1   # a single job for each biological condition

    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='_done_detecting_similar_pssms.txt')

    # 3 steps together!! align each cluster; clean each alignment; calculate pssm for each alignment
    align_clean_pssm_weblogo(biological_conditions, max_msas_per_bc,
                             motif_inference_output_path, logs_dir, error_path, queue_name, verbose)

    # TODO: do the split BEFORE the cutoffs computation so the cutoff computation can be parallelized!!
    # Compute pssm cutoffs for each bc
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: computing pssms cuttoffs for the following biological conditions:\n'
                f'{biological_conditions}')
    script_name = 'calculate_pssm_cutoffs.py'
    num_of_expected_results = 0
    all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        bc_folder = os.path.join(motif_inference_output_path, bc)
        meme_path = os.path.join(bc_folder, 'meme.txt')
        output_path = os.path.join(bc_folder, 'cutoffs.txt')
        done_path = f'{logs_dir}/10_{bc}_done_compute_cutoffs.txt'
        all_cmds_params.append([meme_path, output_path, done_path])

    for cmds_params, bc in zip(all_cmds_params, biological_conditions):
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             [cmds_params],
                             logs_dir, f'{bc}_cutoffs',
                             queue_name, verbose)
        num_of_expected_results += 1  # a single job for each biological condition

    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='_done_compute_cutoffs.txt')

    # Split memes and cutoffs
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: splitting pssms and cuttoffs for paralellizing p-values step:\n'
                f'{biological_conditions}')
    script_name = 'split_meme_and_cutoff_files.py'
    num_of_expected_results = 0
    all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        bc_folder = os.path.join(motif_inference_output_path, bc)
        meme_path = os.path.join(bc_folder, 'meme.txt')
        cutoff_path = os.path.join(bc_folder, 'cutoffs.txt')
        done_path = f'{logs_dir}/11_{bc}_done_split.txt'
        all_cmds_params.append([meme_path, cutoff_path, done_path])

    for cmds_params, bc in zip(all_cmds_params, biological_conditions):
        submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                             [cmds_params],
                             logs_dir, f'{bc}_split',
                             queue_name, verbose)
        num_of_expected_results += 1  # a single job for each biological condition

    wait_for_results(script_name, logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='_done_split.txt')


    # TODO: fix this bug with a GENERAL WRAPPER done_path
    # wait_for_results(script_name, num_of_expected_results)
    with open(motif_inference_done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('parsed_fastq_results', type=str, help='A path in which each subfolder corresponds to a samplename and contains a collapsed faa file')
    parser.add_argument('motif_inference_results', type=str, help='output folder')
    parser.add_argument('logs_dir', type=str, help='logs folder')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='A path to the sample name to biological condition file')
    parser.add_argument('done_file_path', help='A path to a file that signals that the module finished running successfully.')

    parser.add_argument('-q', '--queue', default='pupkoweb', type=str, help='a queue to which the jobs will be submitted')
    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('--max_msas_per_sample', default=100, type=int,
                        help='For each sample, align only the biggest $max_msas_per_sample')
    parser.add_argument('--max_msas_per_bc', default=400, type=int,
                        help='For each biological condition, align only the biggest $max_msas_per_bc')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger('main')

    error_path = args.error_path if args.error_path else os.path.join(args.parsed_fastq_results, 'error.txt')

    infer_motifs(args.parsed_fastq_results, args.max_msas_per_sample, args.max_msas_per_bc,
                 args.motif_inference_results, args.logs_dir,
                 args.samplename2biologicalcondition_path, args.done_file_path,
                 args.queue, True if args.verbose else False, error_path, sys.argv)
