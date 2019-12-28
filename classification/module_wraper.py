import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import *



def build_classifier(first_phase_output_path, motif_inference_output_path,
                     classification_output_path, logs_dir, samplename2biologicalcondition_path,
                     queue_name, verbose, error_path, argv):

    os.makedirs(classification_output_path, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    pvalues_done_path = f'{logs_dir}/pvalues_done.txt'

    samplename2biologicalcondition = load_table(samplename2biologicalcondition_path, 'Barcode {} belongs to more than one sample!!')
    sample_names = sorted(samplename2biologicalcondition)
    biological_conditions = sorted(set(samplename2biologicalcondition.values()))

    # compute pvalues
    logger.info('_'*100)
    logger.info(f'{datetime.datetime.now()}: upper casing all sequences in the faa files')
    script_name = 'upper_case_sequences.py'
    num_of_expected_results = 0
    upper_faa_paths = [] # keep all faas' paths for the next step
    num_of_cmds_per_job = 10
    all_cmds_params = []  # a list of lists. Each sublist contain different parameters set for the same script to reduce the total number of jobs
    for bc in biological_conditions:
        # get current biological condition (splitted) motifs folder
        memes_path = os.path.join(motif_inference_output_path, bc, 'memes')
        for i in range(len(os.listdir(memes_path))):
            # extract each split of the motifs to scan peptides against it
            for sample in sample_names:
                pass #TODO: update from here.
                dir_path = os.path.join(first_phase_output_path, sample_name)
                assert(os.path.exists(dir_path))
                if not os.path.isdir(dir_path):
                    # skip files or folders of non-related biological condition
                    continue
                faa_filename = f'{sample_name}_unique_rpm.faa'
                in_faa_path = os.path.join(first_phase_output_path, sample_name, faa_filename)
                out_faa_dir = os.path.join(classification_output_path, sample_name)
                os.makedirs(out_faa_dir, exist_ok=True)
                out_faa_path = os.path.join(out_faa_dir, f'{sample_name}_upper{faa_filename.split(sample_name)[-1]}')
                upper_faa_paths.append(out_faa_path)
                done_path = f'{logs_dir}/01_{sample_name}_done_uppering.txt'
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
                     error_file_path=error_path, suffix='uppering.txt')


    # TODO: fix this bug with a GENERAL WRAPPER done_path
    # wait_for_results(script_name, num_of_expected_results)
    with open(pvalues_done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('parsed_fastq_results', type=str, help='A path in which each subfolder corresponds to a samplename and contains a collapsed faa file')
    parser.add_argument('motif_inference_results', type=str, help='output folder')
    parser.add_argument('logs_dir', type=str, help='logs folder')
    # parser.add_argument('barcode2samplename', type=str, help='A path to the barcode to sample name file')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='A path to the sample name to biological condition file')
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

    build_classifier(args.parsed_fastq_results, args.max_msas_per_sample, args.max_msas_per_bc,
                 args.motif_inference_results, args.logs_dir,
                 args.samplename2biologicalcondition_path,
                 args.queue, True if args.verbose else False, error_path, sys.argv)
