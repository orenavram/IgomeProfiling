import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import fetch_cmd, wait_for_results
from global_params import src_dir


def run_first_phase(fastq_path, first_phase_output_path, logs_dir, barcode2samplename,
                    left_construct, right_construct, max_mismatches_allowed, min_sequencing_quality,
                    gz, verbose, error_path, argv='no argv'):

    os.makedirs(first_phase_output_path, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    first_phase_done_path = f'{logs_dir}/filter_reads_done.txt'
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
                      f'--min_sequencing_quality {min_sequencing_quality}'] + (['--gz'] if gz else [])

        fetch_cmd(f'{src_dir}/reads_filtration/filter_reads.py',
                  parameters, verbose, error_path)
        num_of_expected_results = 1
        wait_for_results('filter_reads.py', logs_dir, num_of_expected_results,
                         error_file_path=error_path, suffix='demultiplexing.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping filter_reads.py ({done_path} exists)')

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
                if not file.startswith(f'{dir_name}.faa'):  # maybe there's a .gz afterwards
                    continue

                sample_name = file.split('.faa')[0]
                file_path = f'{dir_path}/{file}'
                output_file_path = f'{dir_path}/{sample_name}_unique_rpm.faa'
                done_path = f'{logs_dir}/02_{sample_name}_done_collapsing.txt'

                parameters = [file_path, output_file_path, done_path,
                              '--rpm', f'{first_phase_output_path}/rpm_factors.txt']
                fetch_cmd(f'{src_dir}/reads_filtration/count_and_collapse_duplicates.py',
                          parameters, verbose, error_path)

                # file_path = output_file_path
                # output_file_path = f'{os.path.splitext(file_path)[0]}_cysteine_trimmed.faa'
                # parameters = [file_path, output_file_path]
                # fetch_cmd(f'{src_dir}/reads_filtration/remove_cysteine_loop.py',
                #             parameters, verbose, error_path)
                # TODO: if remove_cysteine_loop.py is fetched, the counts should be recalculated!
                # E.g.: CAAAAC and AAAA are the same after removing Cys

                num_of_expected_results += 1
                break

        wait_for_results('count_and_collapse_duplicates.py', logs_dir, num_of_expected_results,
                     error_file_path=error_path, suffix='collapsing.txt')
        with open(collapsing_done_path, 'w') as f:
            f.write(' '.join(argv) + '\n')


    else:
        logger.info(f'{datetime.datetime.now()}: skipping count_and_collapse_duplicates.py ({done_path} exists)')

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

    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('--left_construct', type=str, default="CAACGTGGC", help='left constant sequence')
    parser.add_argument('--right_construct', type=str, default="GCCT", help='right constant sequence')
    parser.add_argument('--max_mismatches_allowed', type=int, default=1,
                        help='number of mismatches allowed together in both constructs')
    parser.add_argument('--min_sequencing_quality', type=int, default=38,
                        help='Minimum average sequencing threshold allowed after filtration'
                             'for more details, see: https://en.wikipedia.org/wiki/Phred_quality_score')
    parser.add_argument('--gz', action='store_true', help='gzip fastq, filtration_log, fna, and faa files')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    error_path = args.error_path if args.error_path else os.path.join(args.parsed_fastq_results, 'error.txt')

    run_first_phase(args.fastq_path, args.parsed_fastq_results, args.logs_dir,
                    args.barcode2samplename, args.left_construct,
                    args.right_construct, args.max_mismatches_allowed,
                    args.min_sequencing_quality, True if args.gz else False,
                    True if args.verbose else False, error_path, sys.argv)
