import datetime
from operator import ne
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

from auxiliaries.pipeline_auxiliaries import fetch_cmd, wait_for_results, load_table_to_dict
from global_params import src_dir

map_names_json_file = {
    "fastq": "fastq_path",
    "output_reads": "first_phase_output_path",
    "logs_dir": "logs_dir",
    "barcode2sample": "barcode2samplename",
    "done_path": "first_phase_done_path",
    "left_construct": "left_construct",
    "right_construct": "right_construct",
    "max_mismatches_allowed": "max_mismatches_allowed",
    "min_sequencing_quality": "min_sequencing_quality",
    "minimal_length_required": "minimal_length_required",
    "gz": "gz",
    "v": "verbose",
    "error_path": "error_path",
    "queue": "queue",
    "rpm": "rpm",
    "m": "use_mapitope"
}

map_names_command_line = {
    "fastq_path": "fastq_path",
    "parsed_fastq_results": "first_phase_output_path",
    "logs_dir": "logs_dir",
    "barcode2samplename": "barcode2samplename",
    "done_file_path": "first_phase_done_path",
    "left_construct": "left_construct",
    "right_construct": "right_construct",
    "max_mismatches_allowed": "max_mismatches_allowed",
    "min_sequencing_quality": "min_sequencing_quality",
    "minimal_length_required": "minimal_length_required",
    'multi_experiments_config': 'multi_experiments_config',
    "gz": "gz",
    "verbose": "verbose",
    "error_path": "error_path",
    "queue": "queue",
    "rpm": "rpm",
    "mapitope": "use_mapitope"
}

def change_key_name(old_names_dict, map_name):
    new_dict = {}
    if old_names_dict:
        keys = old_names_dict.keys()
        for key in keys:
            new_dict[map_name[key]] = old_names_dict[key]
        return new_dict
    return old_names_dict    


def process_params(args, multi_experiments_config, argv):
    # create data structure for running filter_reads
    base_map =  args.__dict__
    keys = base_map.keys()
    base_map = change_key_name(base_map, map_names_command_line)
    if multi_experiments_config:
        #TODO validation of json file structure    
        f = open(multi_experiments_config)
        multi_experiments_dict = json.load(f)
        configuration = multi_experiments_dict['configuration']
        configuration = change_key_name(configuration, map_names_json_file)
        base_map.update(configuration)
        runs = multi_experiments_dict['runs']
        for run in runs:
            dict_params = base_map.copy()
            exp_params = change_key_name(runs[run], map_names_json_file)
            dict_params.update(exp_params)
            # create new list of argv of the specific run.
            argv_new = []
            argv_new.append(argv[0])
            for k in keys:
                val = str(dict_params[map_names_command_line[k]])
                if (val != 'None') and (val != 'False'):
                    argv_new.append(k)
                    argv_new.append(val)              
            run_first_phase(dict_params['fastq_path'], dict_params['first_phase_output_path'], dict_params['logs_dir'], dict_params['barcode2samplename'], dict_params['first_phase_done_path'],
                    dict_params['left_construct'], dict_params['right_construct'], dict_params['max_mismatches_allowed'], dict_params['min_sequencing_quality'], dict_params['minimal_length_required'],
                    dict_params['rpm'], dict_params['gz'], dict_params['verbose'], dict_params['use_mapitope'], dict_params['error_path'], dict_params['queue'], run, argv_new)    
    else:
        exp_name = ''
        run_first_phase(base_map['fastq_path'], base_map['first_phase_output_path'], base_map['logs_dir'], base_map['barcode2samplename'], base_map['first_phase_done_path'],
                    base_map['left_construct'], base_map['right_construct'], base_map['max_mismatches_allowed'], base_map['min_sequencing_quality'], base_map['minimal_length_required'],
                    base_map['rpm'], base_map['gz'], base_map['verbose'], base_map['use_mapitope'], base_map['error_path'], base_map['queue'], exp_name, argv)


def run_first_phase(fastq_path, first_phase_output_path, logs_dir, barcode2samplename, first_phase_done_path,
                    left_construct, right_construct, max_mismatches_allowed, min_sequencing_quality, minimal_length_required,
                    rpm, gz, verbose, use_mapitope, error_path, queue, exp_name, argv):

    if os.path.exists(first_phase_done_path):
        logger.info(f'{datetime.datetime.now()}: skipping reads_filtration step ({first_phase_done_path} already exists)')
        return        
        
    os.makedirs(first_phase_output_path, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
        
    barcode2samplename_dict = load_table_to_dict(barcode2samplename, 'Barcode {} belongs to more than one sample_name!!')
    sample_names = sorted(barcode2samplename_dict.values())
    error_path = error_path if error_path else os.path.join(first_phase_output_path, 'error.txt')

    script_name = 'filter_reads.py'
    num_of_expected_results = 0
    done_path = f'{logs_dir}/{exp_name}_done_demultiplexing.txt'
    if not os.path.exists(done_path):
        logger.info('_' * 100)
        logger.info(f'{datetime.datetime.now()}: demultiplexig sequences for {done_path}')
        all_cmds_params = []
        if not os.path.exists(done_path):
            cmds = [fastq_path, first_phase_output_path, logs_dir,
                    done_path, barcode2samplename, 'summary_log.txt',
                    f'--error_path {error_path}',
                    f'--left_construct {left_construct}',
                    f'--right_construct {right_construct}',
                    f'--max_mismatches_allowed {max_mismatches_allowed}',
                    f'--min_sequencing_quality {min_sequencing_quality}',
                    f'--minimal_length_required {minimal_length_required}'] + (['--gz'] if gz else [])
            all_cmds_params.append(cmds)
        else:
            logger.debug(f'skipping filter reads as {done_path} found')
            num_of_expected_results += 1
        
        if len(all_cmds_params)>0:
            script_path=f'{src_dir}/reads_filtration/{script_name}'
            for cmds_params in all_cmds_params:
                cmd = fetch_cmd(script_path, cmds_params, verbose, error_path, done_path)
                num_of_expected_results += 1
                
            wait_for_results(script_name, logs_dir, num_of_expected_results,
                                    error_file_path=error_path, suffix=f'_done_demultiplexing.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping filter_reads.py ({done_path} exists)')

    if use_mapitope:
        script_name = 'mapitope_conversion.py'
        mapitope_done_path = f'{logs_dir}/01_done_mapitope_encoding.txt'
        if not os.path.exists(mapitope_done_path):
            logger.info('_' * 100)
            logger.info(f'{datetime.datetime.now()}: mapitope encoding data {mapitope_done_path}')
            # run mapitope_conversion.py
            num_of_expected_results = 0
            for sample_name in sorted(sample_names):
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
            logger.info(f'{datetime.datetime.now()}: skipping mapitope_conversion.py ({mapitope_done_path} exists)')
        
    script_name = 'count_and_collapse_duplicates.py'
    collapsing_done_path = f'{logs_dir}/02_done_collapsing_all.txt'
    if not os.path.exists(collapsing_done_path):
        logger.info('_' * 100)
        logger.info(f'{datetime.datetime.now()}: demultiplexig sequences for {collapsing_done_path}')
        num_of_expected_results = 0
        for dir_name in sorted(sample_names):
            dir_path = os.path.join(first_phase_output_path, dir_name)
            if not os.path.isdir(dir_path):
                continue
            for file in os.listdir(dir_path):
                # look for faa files to collapse
                if not file.startswith(f'{dir_name}.faa') and not file.startswith(f'{dir_name}_mapitope.faa'):  # maybe there's a .gz afterwards
                    continue

                sample_name = file.split('.faa')[0]
                file_path = f'{dir_path}/{file}'
                output_file_path = f'{dir_path}/{sample_name}_unique_rpm.faa'
                done_path = f'{logs_dir}/02_{sample_name}_done_collapsing.txt'
                factors_file_name = 'mapitop_rpm_factors' if 'mapitope' in file else 'rpm_factors'
                
                if not os.path.exists(done_path):
                    if rpm:
                        rpm_factors_path =  f'{first_phase_output_path}/{factors_file_name}.txt'
                        parameters.append('--rpm')
                        parameters.append(rpm_factors_path)        
                        fetch_cmd(f'{src_dir}/reads_filtration/count_and_collapse_duplicates.py', parameters, verbose, error_path, done_path)
                        num_of_expected_results += 1

                else:
                    logger.debug(f'skipping filter reads as {done_path} found')
                    num_of_expected_results += 1    
                

        wait_for_results('count_and_collapse_duplicates.py', logs_dir, num_of_expected_results,
                    error_file_path=error_path, suffix='collapsing.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping count_and_collapse_duplicates.py ({collapsing_done_path} exists)')

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

    parser.add_argument('--multi_experiments_config', type=str, help='A path to json file that contains a match between name of run and fastq, barcode2sample')
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

    process_params(args,args.multi_experiments_config, sys.argv)
