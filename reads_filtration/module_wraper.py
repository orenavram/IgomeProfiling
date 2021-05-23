import datetime
import os
import sys
import json
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import fetch_cmd, wait_for_results, submit_pipeline_step
from global_params import src_dir


def run_first_phase(fastq_path, first_phase_output_path, logs_dir, barcode2samplename, first_phase_done_path,
                    left_construct, right_construct, max_mismatches_allowed, min_sequencing_quality, minimal_length_required,
                    multi_experiments_config, gz, verbose, use_mapitope, error_path, queue, argv='no_argv'):

    if os.path.exists(first_phase_done_path):
        logger.info(f'{datetime.datetime.now()}: skipping reads_filtration step ({first_phase_done_path} already exists)')
        return

    # create data structure for running filter_reads
    map_multi_experiments = {}
    if multi_experiments_config:
        f = open(multi_experiments_config)
        multi_experiments_dict = json.load(f)
        map_multi_experiments = multi_experiments_dict['runs']
    else:
        map_multi_experiments['one'] = {
            "fastq_path": fastq_path,
            "barcode2samplename": barcode2samplename,
            "output_path": first_phase_output_path,
            "logs_dir": logs_dir,
        }
    
    for exp in map_multi_experiments:
        os.makedirs(map_multi_experiments[exp]['output_path'], exist_ok=True)
        os.makedirs(map_multi_experiments[exp]['logs_dir'], exist_ok=True)

    
    script_name = 'filter_reads.py'
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: demultiplexig sequences by {script_name}')
    num_of_expected_results = {} 
    all_cmds_params = []
    for exp in map_multi_experiments:
        map_args_exp = map_multi_experiments[exp]
        done_path=os.path.join(map_args_exp['logs_dir'], f'filter_reads_done.txt')
        if map_args_exp['logs_dir'] not in num_of_expected_results:
            num_of_expected_results[map_args_exp['logs_dir']] = 0
        if not os.path.exists(done_path):
            cmds = [map_args_exp['fastq_path'], map_args_exp['output_path'], map_args_exp['logs_dir'],
                        done_path, map_args_exp['barcode2samplename'], 'summary_log.txt',
                        f'--error_path {error_path}',
                        f'--left_construct {left_construct}',
                        f'--right_construct {right_construct}',
                        f'--max_mismatches_allowed {max_mismatches_allowed}',
                        f'--min_sequencing_quality {min_sequencing_quality}',
                        f'--minimal_length_required {minimal_length_required}'] + (['--gz'] if gz else [])
            all_cmds_params.append(cmds)
        else:
            logger.debug(f'skipping filter reads as {done_path} found')
            num_of_expected_results[map_args_exp['logs_dir']] += 1
    if len(all_cmds_params)>0:
        script_path=f'{src_dir}/reads_filtration/{script_name}'
        for cmds_params, exp in zip(all_cmds_params, map_multi_experiments):
            cmd = fetch_cmd(script_path, cmds_params, verbose, error_path, done_path)
            log_dir = cmds_params[2]
            num_of_expected_results[log_dir] += 1
        wait_for_results(script_name, log_dir, num_of_expected_results[log_dir],
                            error_file_path=error_path, suffix='filter_reads_done.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping filter_reads.py, all reads exists')


    if use_mapitope:
        mapitope_done_path = f'{logs_dir}/01_done_mapitope_encoding.txt'
        if not os.path.exists(mapitope_done_path):
            logger.info('_' * 100)
            logger.info(f'{datetime.datetime.now()}: mapitope encoding data {first_phase_output_path}')
            # run mapitope_conversion.py

            num_of_expected_results = 0
            for sample_name in sorted(os.listdir(first_phase_output_path)):
                dir_path = os.path.join(first_phase_output_path, sample_name)
                if not os.path.isdir(dir_path):
                    continue
                done_path = f'{logs_dir}/01_{sample_name}_done_converting_to_mapitope.txt'
                parameters = [
                    f'{first_phase_output_path}/{sample_name}/{sample_name}.faa',
                    f'{first_phase_output_path}/{sample_name}/{sample_name}_mapitope.faa',
                    done_path
                ]
                fetch_cmd(f'{src_dir}/reads_filtration/mapitope_conversion.py', parameters, verbose, error_path, done_path)
                num_of_expected_results += 1

            wait_for_results('mapitope_conversion.py', logs_dir, num_of_expected_results, error_file_path=error_path, suffix='mapitope.txt')
            with open(mapitope_done_path, 'w') as f:
                f.write(' '.join(argv) + '\n')

        else:
            logger.info(f'{datetime.datetime.now()}: skipping mapitope_conversion.py ({done_path} exists)')
     
    script_name = 'count_and_collapse_duplicates.py'
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: demultiplexig sequences by {script_name}')
    num_of_expected_results = {}
    list_output_path = []
    log_dirs = []
    for exp in map_multi_experiments:
        list_output_path.append(map_multi_experiments[exp]['output_path'])
        log_dirs.append(map_multi_experiments[exp]['logs_dir'])
    list_output_path = list(set(list_output_path))
    log_dirs = list(set(log_dirs))
    for num, output in enumerate(list_output_path):
        if output not in  num_of_expected_results:
            num_of_expected_results[output] = 0
        for dir_name in sorted(os.listdir(output)):
            dir_path = os.path.join(output, dir_name)
            if not os.path.isdir(dir_path):
                continue
            for file in os.listdir(dir_path):
                # look for faa files to collapse
                if not file.startswith(f'{dir_name}.faa') and not file.startswith(f'{dir_name}_mapitope.faa'):  # maybe there's a .gz afterwards
                    continue

                sample_name = file.split('.faa')[0]
                file_path = f'{dir_path}/{file}'
                output_file_path = f'{dir_path}/{sample_name}_unique_rpm.faa'
                done_path = f'{log_dirs[num]}/02_{sample_name}_done_collapsing.txt'
                factors_file_name = 'mapitop_rpm_factors' if 'mapitope' in file else 'rpm_factors'
                rpm_factors_path =  f'{first_phase_output_path}/{factors_file_name}.txt'
                if not os.path.exists(done_path):
                    parameters = [file_path, output_file_path, done_path, '--rpm', rpm_factors_path]
                    fetch_cmd(f'{src_dir}/reads_filtration/count_and_collapse_duplicates.py', parameters, verbose, error_path, done_path)          
                    num_of_expected_results[output] += 1
                else:
                    logger.debug(f'skipping filter reads as {done_path} found')
                    num_of_expected_results[output] += 1    

        wait_for_results('count_and_collapse_duplicates.py', log_dirs[num], num_of_expected_results[output],
                     error_file_path=error_path, suffix='collapsing.txt')


    with open(first_phase_done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}', flush=True)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_path', type=str, help='A fastq file to parse')
    parser.add_argument('parsed_fastq_results', type=str, help='output folder')
    parser.add_argument('logs_dir', type=str, help='logs folder')
    parser.add_argument('barcode2samplename', type=str, help='A path to the barcode to sample name file')
    parser.add_argument('left_construct', default='CAACGTGGC', help='The (constant) sequence from the LEFT of the random sequence') # in exp12: "CAACGTGGC"
    parser.add_argument('right_construct', default='GCCT', help='The (constant) sequence from the RIGHT of the random sequence') # in exp12: "GCCT"
    parser.add_argument('max_mismatches_allowed', type=int, default=1,
                        help='number of mismatches allowed together in both constructs')
    parser.add_argument('min_sequencing_quality', type=int, default=38,
                        help='Minimum average sequencing threshold allowed after filtration'
                             'for more details, see: https://en.wikipedia.org/wiki/Phred_quality_score')
    parser.add_argument('done_file_path', help='A path to a file that signals that the module finished running successfully.')
    parser.add_argument('minimal_length_required', default=3, type=int,
                        help='Shorter peptides will be discarded')

    parser.add_argument('--multi_experiments_config', type=str, help='a path to json file that contains a match between name of run and fastq, barcode2sample')
    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('--gz', action='store_true', help='gzip fastq, filtration_log, fna, and faa files')
    parser.add_argument('-q', '--queue', default='pupkoweb', type=str, help='a queue to which the jobs will be submitted')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    parser.add_argument('-m', '--mapitope', action='store_true', help='use mapitope encoding')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    error_path = args.error_path if args.error_path else os.path.join(args.parsed_fastq_results, 'error.txt')

    run_first_phase(args.fastq_path, args.parsed_fastq_results, args.logs_dir,
                    args.barcode2samplename, args.done_file_path, args.left_construct,
                    args.right_construct, args.max_mismatches_allowed,
                    args.min_sequencing_quality, args.minimal_length_required, args.multi_experiments_config, True if args.gz else False,
                    True if args.verbose else False, args.mapitope, error_path, args.queue, sys.argv)
