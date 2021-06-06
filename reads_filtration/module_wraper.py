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

from auxiliaries.pipeline_auxiliaries import fetch_cmd, wait_for_results
from auxiliaries.stop_machine_aws import stop_machines
from auxiliaries.validation_files import is_input_files_valid 

def run_first_phase(fastq_path, first_phase_output_path, logs_dir, barcode2samplename, first_phase_done_path,
                    left_construct, right_construct, max_mismatches_allowed, min_sequencing_quality, minimal_length_required,
                    check_files_valid, stop_machines_flag, type_machines_to_stop, name_machines_to_stop,
                    rpm, gz, verbose, use_mapitope, error_path, queue, argv='no_argv'):
  
    # check the validation of files barcode2samplename_path and samplename2biologicalcondition_path
    if check_files_valid and not is_input_files_valid(samplename2biologicalcondition_path='', barcode2samplename_path=barcode2samplename, logger=logger):
        return

    os.makedirs(first_phase_output_path, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    if os.path.exists(first_phase_done_path):
        logger.info(f'{datetime.datetime.now()}: skipping reads_filtration step ({first_phase_done_path} already exists)')
        return

    done_path = f'{logs_dir}/done_demultiplexing.txt'
    logger.info('_' * 100)
    logger.info(f'{datetime.datetime.now()}: demultiplexig sequences for {first_phase_output_path}')
    if not os.path.exists(done_path):
        # run filter_reads.py
        parameters = [fastq_path, first_phase_output_path, logs_dir,
                      done_path, barcode2samplename,
                      f'--error_path {error_path}',
                      f'--left_construct {left_construct}',
                      f'--right_construct {right_construct}',
                      f'--max_mismatches_allowed {max_mismatches_allowed}',
                      f'--min_sequencing_quality {min_sequencing_quality}',
                      f'--minimal_length_required {minimal_length_required}'] + (['--gz'] if gz else [])

        fetch_cmd(f'{src_dir}/reads_filtration/filter_reads.py',
                  parameters, verbose, error_path)
        num_of_expected_results = 1
        wait_for_results('filter_reads.py', logs_dir, num_of_expected_results,
                         error_file_path=error_path, suffix='demultiplexing.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping filter_reads.py ({done_path} exists)')

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

    collapsing_done_path = f'{logs_dir}/02_done_collapsing_all.txt'
    if not os.path.exists(collapsing_done_path):
        logger.info('_' * 100)
        logger.info(f'{datetime.datetime.now()}: counting and collapsing duplicated sequences for {first_phase_output_path}')
        # run count_and_collapse_duplicates.py and remove_cysteine_loop.py
        num_of_expected_results = 0
        for dir_name in sorted(os.listdir(first_phase_output_path)):
            dir_path = os.path.join(first_phase_output_path, dir_name)
            if not os.path.isdir(dir_path):
                continue
            for file in os.listdir(dir_path):
                # look for faa files to collapse
                if not file.startswith(f'{dir_name}.faa') and not file.startswith(f'{dir_name}_mapitope.faa'): # maybe there's a .gz afterwards
                    continue

                sample_name = file.split('.faa')[0]
                file_path = f'{dir_path}/{file}'
                output_file_path = f'{dir_path}/{sample_name}_unique_rpm.faa'
                done_path = f'{logs_dir}/02_{sample_name}_done_collapsing.txt'
                factors_file_name = 'mapitop_rpm_factors' if 'mapitope' in file else 'rpm_factors'
                parameters = [file_path, output_file_path, done_path]
                if rpm:
                    rpm_factors_path =  f'{first_phase_output_path}/{factors_file_name}.txt'
                    parameters.append('--rpm')
                    parameters.append(rpm_factors_path)

                fetch_cmd(f'{src_dir}/reads_filtration/count_and_collapse_duplicates.py', parameters, verbose, error_path, done_path)

                num_of_expected_results += 1
                break

        wait_for_results('count_and_collapse_duplicates.py', logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='collapsing.txt')
        with open(collapsing_done_path, 'w') as f:
            f.write(' '.join(argv) + '\n')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping count_and_collapse_duplicates.py ({done_path} exists)')

    if stop_machines_flag:
        stop_machines(type_machines_to_stop, name_machines_to_stop, logger)
        
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
    parser.add_argument('done_file_path', help='A path to a file that signals that the module finished running successfully')
    parser.add_argument('minimal_length_required', default=3, type=int,
                        help='Shorter peptides will be discarded')
    
    parser.add_argument('--check_files_valid', action='store_true', help='Need to check the validation of the files (samplename2biologicalcondition_path / barcode2samplenaem)')
    parser.add_argument('--stop_machines', action='store_true', help='Turn off the machines in AWS at the end of the running')
    parser.add_argument('--type_machines_to_stop', default='', type=str, help='Type of machines to stop, separated by comma. Empty value means all machines. Example: t2.2xlarge,m5a.24xlarge')
    parser.add_argument('--name_machines_to_stop', default='', type=str, help='Names (patterns) of machines to stop, separated by comma. Empty value means all machines. Example: worker*')
    parser.add_argument('--rpm', action='store_true', help='Normalize counts to "reads per million" (sequence proportion x 1,000,000)')
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
                    args.min_sequencing_quality, args.minimal_length_required,
                    args.check_files_valid, args.stop_machines, args.type_machines_to_stop, args.name_machines_to_stop,
                    args.rpm, args.gz, args.verbose, args.mapitope, error_path, args.queue, sys.argv)
