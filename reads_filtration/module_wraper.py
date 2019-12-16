import os
from Auxiliaries.pipeline_auxiliaries import fetch_cmd
from global_params import src_dir


def run_first_phase(fastq_path, parsed_fastq_results, barcode2samplename,
                    left_construct, right_construct, max_mismatches_allowed, min_sequencing_quality,
                    gz, verbose):

    results_output = f'{parsed_fastq_results}/results/'
    # logs_output = f'{parsed_fastq_results}/logs/'
    done_path = f'{results_output}/done.txt'

    if not os.path.exists(done_path):
        # run filter_reads.py
        parameters = [fastq_path, results_output, barcode2samplename,
                      f'--left_construct {left_construct}',
                      f'--right_construct {right_construct}',
                      f'--max_mismatches_allowed {max_mismatches_allowed}',
                      f'--min_sequencing_quality {min_sequencing_quality}']
        if gz:
            parameters.append('--gz')
        if verbose:
            parameters.append('-v')

        fetch_cmd(f'{src_dir}/reads_filtration/filter_reads.py', parameters)

    # run count_and_collapse_duplicates.py and remove_cysteine_loop.py
    for dir_name in sorted(os.listdir(results_output)):
        dir_path = os.path.join(results_output, dir_name)
        if not os.path.isdir(dir_path):
            continue
        for file in os.listdir(dir_path):
            # look for faa files to collapse
            if not file.startswith(f'{dir_name}.faa'):  # maybe there's a .gz afterwards
                continue

            sample_name = file.split('.faa')[0]
            file_path = f'{dir_path}/{file}'
            output_file_path = f'{dir_path}/{sample_name}_unique_rpm.faa'
            parameters = [file_path, output_file_path, '--rpm', f'{results_output}/rpm_factors.txt']
            fetch_cmd(f'{src_dir}/reads_filtration/count_and_collapse_duplicates.py', parameters)

            # file_path = output_file_path
            # output_file_path = f'{os.path.splitext(file_path)[0]}_cysteine_trimmed.faa'
            # parameters = [file_path, output_file_path]
            # fetch_cmd(f'{src_dir}/reads_filtration/remove_cysteine_loop.py', parameters)
            # TODO: if remove_cysteine_loop.py is fetched, the counts should be recalculated!
            # E.g.: CAAAAC and AAAA are the same after removing Cys

            # num_of_expected_results += 1
            break

    # wait_for_results(script_name, num_of_expected_results)
    with open(done_path, 'w'):
        pass




if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_path', type=str, help='A fastq file to parse')
    parser.add_argument('parsed_fastq_results', type=str, help='folder output')
    parser.add_argument('barcode2samplename', type=str, help='A path to the barcode to sample name file')

    parser.add_argument('--left_construct', type=str, default="CAACGTGGC", help='SfilSite')
    parser.add_argument('--right_construct', type=str, default="GCCT", help='RightConstruct')
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

    run_first_phase(args.fastq_path, args.parsed_fastq_results,
                    args.barcode2samplename, args.left_construct,
                    args.right_construct, args.max_mismatches_allowed,
                    args.min_sequencing_quality, True if args.gz else False,
                    True if args.verbose else False)
