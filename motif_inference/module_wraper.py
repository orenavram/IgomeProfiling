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


def align_clean_pssm_weblogo(folder_names_to_handle, max_clusters_to_align, gap_frequency,
                             motif_inference_output_path, logs_dir, error_path, queue_name, verbose, data_type):
    # For each sample, align each cluster
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: aligning clusters in each sample')
    script_name = 'align_sequences.py'
    done_path_suffix = f'done_msa_{data_type}.txt'
    num_of_expected_results = 0
    msas_paths = []  # keep all msas' paths for the next step
    num_of_cmds_per_job = 33
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    max_number_of_leading_zeros = len(str(max_clusters_to_align))
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
            aligned_cluster_path = os.path.join(aligned_sequences_path, faa_filename)
            done_path = f'{logs_dir}/05_{sample_name}_{cluster_rank}_{done_path_suffix}'
            if not os.path.exists(done_path):
                all_cmds_params.append([unaligned_cluster_path, aligned_cluster_path, done_path])
            else:
                logger.debug(f'Skipping sequence alignment as {done_path} exists')
                num_of_expected_results += 1

    if len(all_cmds_params) > 0:
        for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
            current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
            sample_name = current_batch[0][1].split('/')[-3]
            assert sample_name in folder_names_to_handle, f'Sample {sample_name} not in folder names list:\n{folder_names_to_handle}'
            cluster_rank = get_cluster_rank_from(current_batch[-1][1])
            cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                current_batch,
                                logs_dir, f'{sample_name}_{cluster_rank}_msa',
                                queue_name, verbose)

            num_of_expected_results += len(current_batch)


        wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path=error_path, suffix=done_path_suffix)
    else:
        logger.info('Skipping sequence alignment, all found')


    # For each sample, clean alignments from gappy columns
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: cleaning alignments from gappy columns')
    script_name = 'remove_gappy_columns.py'
    done_path_suffix = f'done_cleaning_msa_{data_type}.txt'
    num_of_expected_results = 0
    cleaned_msas_paths = []  # keep all cleaned msas' paths for the next step
    num_of_cmds_per_job = 50  # a super fast script. No point to put less than 50 (as the overhead will take longer)..
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
            done_path = f'{logs_dir}/06_{sample_name}_{msa_name}_{done_path_suffix}'
            if not os.path.exists(done_path):
                all_cmds_params.append([msa_path, cleaned_msa_path, done_path,
                                        '--maximal_gap_frequency_allowed_per_column', gap_frequency])
            else:
                logger.debug(f'Skipping cleaning as {done_path} exists')
                num_of_expected_results += 1

    if len(all_cmds_params) > 0:
        for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
            current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
            sample_name = current_batch[0][1].split('/')[-3]
            assert sample_name in folder_names_to_handle, f'Sample {sample_name} not in folder names list:\n{folder_names_to_handle}'
            cluster_rank = get_cluster_rank_from(current_batch[-1][0])
            cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                current_batch,
                                logs_dir, f'{sample_name}_{cluster_rank}_clean',
                                queue_name, verbose)
            num_of_expected_results += len(current_batch)

        wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path=error_path, suffix=done_path_suffix) #, done_files_list=done_files_list)
    else:
        logger.info(f'Skipping cleaning, all exists')


    # For each sample, generate a meme file with a corresponding pssm for each alignment
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: generating meme files for each sample from cleaned alignments')
    script_name = 'create_meme.py'
    meme_done_path_suffix = f'done_meme_{data_type}.txt'
    num_of_expected_memes = 0
    num_of_cmds_per_job = 1
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for cleaned_msas_path in cleaned_msas_paths:
        sample_motifs_dir = os.path.split(cleaned_msas_path)[0]
        sample_name = os.path.split(sample_motifs_dir)[-1]
        assert sample_name in folder_names_to_handle, f'Sample {sample_name} not in folder names list:\n{folder_names_to_handle}'
        meme_path = os.path.join(sample_motifs_dir, 'meme.txt')
        done_path = f'{logs_dir}/07_{sample_name}_{meme_done_path_suffix}'
        if not os.path.exists(done_path):
            all_cmds_params.append([cleaned_msas_path, meme_path, done_path])
        else:
            logger.debug(f'Skipping meme creation as {done_path} exists')
            num_of_expected_results += 1
    
    if len(all_cmds_params) > 0:
        for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
            current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
            memes_cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
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
            # all_cmds_params.append([msa_path, weblogo_prefix_path])

    if len(all_cmds_params) > 0:
        for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
            current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
            # submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
            #                      current_batch,
            #                      logs_dir, f'weblogo_{i}th_batch',
            #                      queue_name='pupkolab', verbose=verbose)

    # wait for the memes!! (previous logical block!)
    # (no need to wait for weblogos..)
    if num_of_expected_memes > 0:
        script_name = 'create_meme.py'
        wait_for_results(script_name, logs_dir, num_of_expected_memes, example_cmd=memes_cmd,
                        error_file_path=error_path, suffix=meme_done_path_suffix)
    else:
        logger.info('Skipping memes creation, all found')


def compute_cutoffs_then_split(biological_conditions, meme_split_size,
    motif_inference_output_path, logs_dir, queue_name, verbose):
    # Compute pssm cutoffs for each bc
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: computing pssms cutoffs for the following biological conditions:\n'
                f'{biological_conditions}')
    script_name = 'calculate_pssm_cutoffs.py'
    num_of_expected_results = 0
    all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        bc_folder = os.path.join(motif_inference_output_path, bc)
        meme_path = os.path.join(bc_folder, 'meme.txt')
        output_path = os.path.join(bc_folder, 'cutoffs.txt')
        done_path = f'{logs_dir}/13_{bc}_done_compute_cutoffs.txt'
        if not os.path.exists(done_path):
            all_cmds_params.append([meme_path, output_path, done_path, '--total_memes', 0])
        else:
            logger.debug(f'Skipping cutoff as {done_path} exists')
            num_of_expected_results += 1

    if len(all_cmds_params) > 0:
        for cmds_params, bc in zip(all_cmds_params, biological_conditions):
            cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                [cmds_params],
                                logs_dir, f'{bc}_cutoffs',
                                queue_name, verbose)
            num_of_expected_results += 1  # a single job for each biological condition

        wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path=error_path, suffix='_done_compute_cutoffs.txt')
    else:
        logger.info('Skipping cuttoffs, all exist')

    # Split memes and cutoffs
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: splitting pssms and cutoffs for paralellizing p-values step:\n'
                f'{biological_conditions}')
    script_name = 'split_meme_and_cutoff_files.py'
    num_of_expected_results = 0
    all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        bc_folder = os.path.join(motif_inference_output_path, bc)
        meme_path = os.path.join(bc_folder, 'meme.txt')
        cutoff_path = os.path.join(bc_folder, 'cutoffs.txt')
        done_path = f'{logs_dir}/14_{bc}_done_split.txt'
        if not os.path.exists(done_path):
            all_cmds_params.append([meme_path, cutoff_path, done_path, '--motifs_per_file', meme_split_size])
        else:
            logger.debug(f'Skipping split as {done_path} exists')
            num_of_expected_results += 1

    if len(all_cmds_params) > 0:
        for cmds_params, bc in zip(all_cmds_params, biological_conditions):
            cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                [cmds_params],
                                logs_dir, f'{bc}_split',
                                queue_name, verbose)
            num_of_expected_results += 1  # a single job for each biological condition

        wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path=error_path, suffix='_done_split.txt')
    else:
        logger.info('Skipping split, all exist')


def split_then_compute_cutoffs(biological_conditions, meme_split_size,
    motif_inference_output_path, logs_dir, queue_name, verbose):
    # Count memes per BC (synchrnous)
    memes_per_bc = {}
    for bc in biological_conditions:
        bc_folder = os.path.join(motif_inference_output_path, bc)
        meme_path = os.path.join(bc_folder, 'meme.txt')
        memes_count = count_memes(meme_path)
        memes_per_bc[bc] = memes_count

    # Split memes and cutoffs
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: splitting pssms and cutoffs for paralellizing p-values step:\n'
                f'{biological_conditions}')
    script_name = 'split_meme_and_cutoff_files.py'
    num_of_expected_results = 0
    all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        bc_folder = os.path.join(motif_inference_output_path, bc)
        meme_path = os.path.join(bc_folder, 'meme.txt')
        cutoff_path = 'skip'
        done_path = f'{logs_dir}/13_{bc}_done_split.txt'
        if not os.path.exists(done_path):
            all_cmds_params.append([meme_path, cutoff_path, done_path, '--motifs_per_file', meme_split_size])
        else:
            logger.debug(f'Skipping split as {done_path} exists')
            num_of_expected_results += 1

    if len(all_cmds_params) > 0:
        for cmds_params, bc in zip(all_cmds_params, biological_conditions):
            cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                [cmds_params],
                                logs_dir, f'{bc}_split',
                                queue_name, verbose)
            num_of_expected_results += 1  # a single job for each biological condition

        wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path=error_path, suffix='_done_split.txt')
    else:
        logger.info('Skipping split, all exist')
    
    # Compute pssm cutoffs for each bc
    # TODO change to read to cut files (read directory)
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: computing pssms cutoffs for the following biological conditions:\n'
                f'{biological_conditions}')
    script_name = 'calculate_pssm_cutoffs.py'
    num_of_expected_results = 0
    all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    cutoffs_bcs = []
    for bc in biological_conditions:
        bc_folder = os.path.join(motif_inference_output_path, bc)
        bc_memes_folder = os.path.join(bc_folder, 'memes')
        bc_cutoffs_folder = os.path.join(bc_folder, 'cutoffs')
        os.makedirs(bc_cutoffs_folder, exist_ok=True)
        for file_name in sorted(os.listdir(bc_memes_folder)):
            meme_path = os.path.join(bc_memes_folder, file_name)
            output_path = os.path.join(bc_cutoffs_folder, file_name)
            done_path = f'{logs_dir}/14_{bc}_{file_name}_done_compute_cutoffs.txt'
            if not os.path.exists(done_path):
                all_cmds_params.append([meme_path, output_path, done_path, '--total_memes', memes_per_bc[bc]])
                cutoffs_bcs.append(bc)
            else:
                logger.debug(f'Skipping calculate cutoff as {done_path} exists')
                num_of_expected_results += 1

    if len(all_cmds_params) > 0:
        for cmds_params, bc in zip(all_cmds_params, cutoffs_bcs):
            cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                [cmds_params],
                                logs_dir, f'{bc}_cutoffs',
                                queue_name, verbose)
            num_of_expected_results += 1  # a single job for each biological condition

        wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path=error_path, suffix='_done_compute_cutoffs.txt')
    else:
        logger.info('Skipping calculate cutoffs, all exists')


def infer_motifs(first_phase_output_path, max_msas_per_sample, max_msas_per_bc,
                 max_number_of_cluster_members_per_sample, max_number_of_cluster_members_per_bc,
                 gap_frequency, motif_inference_output_path, logs_dir, samplename2biologicalcondition_path,
                 motif_inference_done_path, multi_experiments_config, threshold, word_length, discard, queue_name, verbose, concurrent_cutoffs, meme_split_size, error_path, use_mapitope, argv):

    
    # create data structure for running infer_motifs
    map_multi_experiments = {}
    if multi_experiments_config:
        f = open(multi_experiments_config)
        multi_experiments_dict = json.load(f)
        map_multi_experiments = multi_experiments_dict['runs']
    else:
        map_multi_experiments['main'] = {
            "samplename2biologicalcondition_path": samplename2biologicalcondition_path,
            "output_path_phase1": first_phase_output_path,
            "output_path_phase2": motif_inference_output_path,
            "logs_dir": logs_dir,
            "second_phase_done_path":motif_inference_done_path
        }
    
    for exp in map_multi_experiments:
        second_phase_done_path_exp = map_multi_experiments[exp]['second_phase_done_path']
        if os.path.exists(second_phase_done_path_exp):
            logger.info(f'{datetime.datetime.now()}: skipping motif_inderence step ({second_phase_done_path_exp} already exists)')
            del map_multi_experiments[exp]
        
        os.makedirs(map_multi_experiments[exp]['output_path_phase2'], exist_ok=True)
        os.makedirs(map_multi_experiments[exp]['logs_dir'], exist_ok=True)
    
    if not map_multi_experiments:
        return

    for exp in map_multi_experiments:
        logger.info(f'{datetime.datetime.now()}: start motif_inference step for experiment: ({exp})')
        exp_params = map_multi_experiments[exp]
        samplename2biologicalcondition = load_table_to_dict(exp_params['samplename2biologicalcondition_path'], 'Barcode {} belongs to more than one sample!!')
        sample_names = sorted(samplename2biologicalcondition)
        biological_conditions = sorted(set(samplename2biologicalcondition.values()))
        logs_dir = exp_params['logs_dir']
        first_phase_output = exp_params['output_path_phase1']
        motif_inference_output_path = exp_params['output_path_phase2']
        # Make sure all sequences are in upper case letters (for example, no need to differentiate between Q and q)
        logger.info('_'*100)
        logger.info(f'{datetime.datetime.now()}: upper casing all sequences in the faa files')
        script_name = 'upper_case_sequences.py'
        num_of_expected_results = 0
        upper_faa_paths = [] # keep all faas' paths for the next step
        for sample_name in sample_names:
            dir_path = os.path.join(first_phase_output, sample_name)
            assert os.path.exists(dir_path), f'reads filtration directory does not exist!\n{dir_path}'
            if not os.path.isdir(dir_path):
                # skip files or folders of non-related biological condition
                continue

            for file_name in os.listdir(dir_path):

                if file_name.endswith('faa') and 'unique' in file_name and ('mapitope' in file_name) == use_mapitope:
                    faa_filename = file_name
                    break
            else:
                raise ValueError(f'No faa file at {dir_path}')

            in_faa_path = os.path.join(first_phase_output, sample_name, faa_filename)
            out_faa_dir = os.path.join(motif_inference_output_path, sample_name)
            os.makedirs(out_faa_dir, exist_ok=True)
            out_faa_path = os.path.join(out_faa_dir, f'{sample_name}_upper{faa_filename.split("unique")[-1]}') # not unique anymore (q->Q)
            upper_faa_paths.append(out_faa_path)
            logs_dir = exp_params['logs_dir']
            done_path = f'{logs_dir}/01_{sample_name}_done_uppering.txt'
            fetch_cmd(f'{src_dir}/motif_inference/{script_name}',
                    [in_faa_path, out_faa_path, done_path],
                    verbose, error_path, done_path)
            num_of_expected_results += 1

        wait_for_results(script_name, logs_dir, num_of_expected_results,
                        error_file_path=error_path, suffix='uppering.txt')

        # Remove flanking Cysteines before clustering
        logger.info('_'*100)
        logger.info(f'{datetime.datetime.now()}: removing flanking Cysteines from faa files')
        script_name = 'remove_cysteine_loop.py'
        num_of_expected_results = 0
        no_cys_faa_paths = []  # keep all faas' paths for the next step
        for upper_faa_path in upper_faa_paths:
            no_cys_faa_path = upper_faa_path.replace('_upper', '_upper_cysteineless')
            no_cys_faa_paths.append(no_cys_faa_path)
            # ~/igomeProfilingPipeline/experiments/exp12/analysis/motif_inference/17b_01/17b_01_upper_unique_rpm.faa
            sample_name = upper_faa_path.split('/')[-1].split('_upper_')[0]
            done_path = f'{logs_dir}/02_{sample_name}_remove_cysteines.txt'
            fetch_cmd(f'{src_dir}/motif_inference/{script_name}',
                    [upper_faa_path, no_cys_faa_path, done_path],
                    verbose, error_path, done_path)
            num_of_expected_results += 1

        wait_for_results(script_name, logs_dir, num_of_expected_results,
                        error_file_path=error_path, suffix='remove_cysteines.txt')

        # Clustering sequences within each sample
        logger.info('_'*100)
        logger.info(f'{datetime.datetime.now()}: clustering sequences in each sample')
        script_name = 'cluster.py'
        num_of_expected_results = 0
        clstr_paths = [] # keep all clusters' paths for the next step
        num_of_cmds_per_job = 1
        all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for no_cys_faa_path in no_cys_faa_paths:
            faa_dir, faa_filename = os.path.split(no_cys_faa_path)
            sample_name = os.path.split(faa_dir)[-1]
            assert sample_name in sample_names, f'Sample {sample_name} not in sample names list:\n{sample_names}'
            output_prefix = os.path.join(faa_dir, sample_name)
            clstr_paths.append(f'{output_prefix}.clstr')
            done_path = f'{logs_dir}/03_{sample_name}_done_clustering.txt'
            if not os.path.exists(done_path):
                cmds = [no_cys_faa_path, output_prefix, done_path , '--threshold', threshold, '--word_length', word_length, '--discard', discard]
                all_cmds_params.append(cmds)
            else:
                logger.debug(f'Skipping clustering as {done_path} exists')
                num_of_expected_results += 1

        if len(all_cmds_params) > 0:
            for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
                current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
                sample_name = os.path.split(current_batch[0][1])[-1]
                assert sample_name in sample_names, f'Sample {sample_name} not in sample names list:\n{sample_names}'
                cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                    current_batch,
                                    logs_dir, f'{sample_name}_cluster', queue_name, verbose)
                num_of_expected_results += len(current_batch)

            wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                            error_file_path=error_path, suffix='clustering.txt')
        else:
            logger.info('Skipping clustering, all exists')

        # For each sample, split the faa file to the clusters inferred in the previous step
        # this step uses the sequences WITH THE FLANKING CYSTEINE so the msa will use these Cs
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
            done_path = f'{logs_dir}/04_{sample_name}_done_extracting_sequences.txt'
            if not os.path.exists(done_path):
                all_cmds_params.append([upper_faa_path, clstr_path, clusters_sequences_path, done_path,
                                        f'--max_num_of_sequences_to_keep {max_number_of_cluster_members_per_sample}',
                                        f'--file_prefix {sample_name}'])
            else:
                logger.debug(f'Skipping sequence extraction as {done_path} exists')
                num_of_expected_results += 1

        if len(all_cmds_params) > 0:
            for i in range(0, len(all_cmds_params), num_of_cmds_per_job):
                current_batch = all_cmds_params[i: i + num_of_cmds_per_job]
                clusters_sequences_path = current_batch[0][2]
                sample_name = clusters_sequences_path.split('/')[-2]
                assert sample_name in sample_names, f'Sample {sample_name} not in sample names list:\n{sample_names}'
                cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                    current_batch,
                                    logs_dir, f'{sample_name}_extracting_sequences', queue_name, verbose)
                num_of_expected_results += len(current_batch)

            wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                            error_file_path=error_path, suffix='extracting_sequences.txt')
        else:
            logger.info('Skipping sequences extraction, all exist')

        # 3 steps together!! align each cluster; clean each alignment; calculate pssm for each alignment
        align_clean_pssm_weblogo(sample_names, max_msas_per_sample, gap_frequency,
                                motif_inference_output_path, logs_dir, error_path, queue_name, verbose, 'samples')

        # Merge memes of the same biological condition
        logger.info('_'*100)
        logger.info(f'{datetime.datetime.now()}: merging meme files for each of the following biological conditions\n'
                    f'{biological_conditions}')
        script_name = 'merge_meme_files.py'
        num_of_expected_results = 0
        biological_condition_memes = []
        all_cmds_params = [] # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
        for bc in biological_conditions:
            done_path = f'{logs_dir}/08_{bc}_done_meme_merge.txt'
            bc_folder = os.path.join(motif_inference_output_path, bc)
            os.makedirs(bc_folder, exist_ok=True)
            output_path = os.path.join(bc_folder, 'merged_meme_sorted.txt')
            biological_condition_memes.append(output_path)
            if not os.path.exists(done_path):
                all_cmds_params.append([motif_inference_output_path, bc, output_path, done_path])
            else:
                logger.debug(f'Skipping merge as {done_path} exists')
                num_of_expected_results += 1

        if len(all_cmds_params) > 0:
            for cmds_params, bc in zip(all_cmds_params, biological_conditions):
                cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                    [cmds_params],
                                    logs_dir, f'{bc}_merge_meme',
                                    queue_name, verbose)
                num_of_expected_results += 1  # a single job for each biological condition

            wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                            error_file_path=error_path, suffix='_done_meme_merge.txt')
        else:
            logger.info('Skipping merge, all exist')


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
            done_path = f'{logs_dir}/09_{bc}_done_detecting_similar_pssms.txt'
            if not os.path.exists(done_path):
                all_cmds_params.append([motif_inference_output_path, merged_meme_path, bc, relevant_samples,
                                        max_number_of_cluster_members_per_bc,
                                        output_path, done_path])
            else:
                logger.debug(f'Skipping unite as {done_path} exists')
                num_of_expected_results += 1

        if len(all_cmds_params) > 0:
            for cmds_params, bc in zip(all_cmds_params, biological_conditions):
                cmd = submit_pipeline_step(f'{src_dir}/motif_inference/{script_name}',
                                    [cmds_params],
                                    logs_dir, f'{bc}_detect_similar_pssms',
                                    queue_name, verbose)
                num_of_expected_results += 1   # a single job for each biological condition

            wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                            error_file_path=error_path, suffix='_done_detecting_similar_pssms.txt')
        else:
            logger.info('Skipping unite, all exist')

        # 3 steps together!! align each cluster; clean each alignment; calculate pssm for each alignment
        align_clean_pssm_weblogo(biological_conditions, max_msas_per_bc, gap_frequency,
                                motif_inference_output_path, logs_dir, error_path, queue_name, verbose, 'biological_conditions')

        if concurrent_cutoffs:
            split_then_compute_cutoffs(biological_conditions, meme_split_size,
                motif_inference_output_path, logs_dir, queue_name, verbose)
        else:
            compute_cutoffs_then_split(biological_conditions, meme_split_size,
                motif_inference_output_path, logs_dir, queue_name, verbose)

        # TODO: fix this bug with a GENERAL WRAPPER done_path
        # wait_for_results(script_name, num_of_expected_results)
        motif_inference_done_path = map_multi_experiments[exp]['second_phase_done_path']
        if exp != 'main':
            argv = ''
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
    parser.add_argument('max_msas_per_sample', type=int,
                        help='For each sample, align only the biggest $max_msas_per_sample')
    parser.add_argument('max_msas_per_bc', type=int,
                        help='For each biological condition, align only the biggest $max_msas_per_bc')
    parser.add_argument('max_number_of_cluster_members_per_sample', type=int,
                        help='How many members (at most) should be taken to each cluster in each sample')
    parser.add_argument('max_number_of_cluster_members_per_bc', type=int,
                        help='How many members (at most) should be taken to each cluster after motif unification')
    parser.add_argument('allowed_gap_frequency',
                        help='Maximal gap frequency allowed in msa (higher frequency columns are removed)',
                        type=lambda x: float(x) if 0 < float(x) < 1
                                                else parser.error(f'The threshold of the maximal gap frequency allowed per column should be between 0 to 1'))

    parser.add_argument('done_file_path', help='A path to a file that signals that the module finished running successfully.')
    
    parser.add_argument('--multi_experiments_config', type=str, help='A path to json file that contains a match between name of run to params to phase 2')
    parser.add_argument('--threshold', default='0.5', help='Minimal sequence similarity threshold required',
                        type=lambda x: float(x) if 0.4 <= float(x) <= 1
                                                else parser.error(f'CD-hit allows thresholds between 0.4 to 1'))
    parser.add_argument('--word_length', default='2', choices=['2', '3', '4', '5'],
                        help='A heuristic of CD-hit. Choose of word size:\n5 for similarity thresholds 0.7 ~ 1.0\n4 for similarity thresholds 0.6 ~ 0.7\n3 for similarity thresholds 0.5 ~ 0.6\n2 for similarity thresholds 0.4 ~ 0.5')
    parser.add_argument('--discard', default='1', help='Include only sequences longer than <$discard> for the analysis. (CD-hit uses only sequences that are longer than 10 amino acids. When the analysis includes shorter sequences, this threshold should be lowered. Thus, it is set here to 1 by default.)')
    parser.add_argument('--concurrent_cutoffs', action='store_true',
                        help='Use new method which splits meme before cutoffs and runs cutoffs concurrently')
    parser.add_argument('--meme_split_size', type=int, default=5,
                        help='Split size, how many meme per files for calculations')
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
    concurrent_cutoffs = True if args.concurrent_cutoffs else False

    infer_motifs(args.parsed_fastq_results, args.max_msas_per_sample, args.max_msas_per_bc,
                 args.max_number_of_cluster_members_per_sample, args.max_number_of_cluster_members_per_bc,
                 args.allowed_gap_frequency, args.motif_inference_results, args.logs_dir, args.samplename2biologicalcondition_path,
                 args.done_file_path, args.multi_experiments_config, args.threshold, args.word_length, args.discard, args.queue,
                 args.verbose, concurrent_cutoffs, args.meme_split_size, error_path, args.mapitope, sys.argv)
