import datetime
import gzip
import os
import sys
import regex
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

# needs $src_dir in path
from auxiliaries.pipeline_auxiliaries import load_table_to_dict, fail


def get_barcodes_dictionaries(barcode_to_samplename, output_dir, gz) -> {str: {str: str}}:
    barcode2filepaths = {}
    barcode2filehandlers = {}
    barcode2statistics = {}
    barcode2info_filenames = {}

    open_function = open
    mode = 'w'
    extension = ''
    if gz:
        open_function = gzip.open
        mode = 'wt'
        extension = '.gz'

    for barcode in barcode_to_samplename:
        sample_name = barcode_to_samplename[barcode]  # f'{barcode}_{os.path.split(fastq_file)[-1]}'

        # create a sub dir for that barcode (under output_dir)
        subdir_path = f'{output_dir}/{sample_name}'
        os.makedirs(subdir_path, exist_ok=True)

        # add filenames
        barcode2filepaths[barcode] = {}
        barcode2filepaths[barcode]['filtration_log'] = f'{subdir_path}/{sample_name}.filtration_log.txt{extension}'

        # add filehandlers
        barcode2filehandlers[barcode] = {}
        barcode2filehandlers[barcode]['filtration_log'] = open_function(barcode2filepaths[barcode]['filtration_log'], mode)

        barcode2statistics[barcode] = dict.fromkeys(['legal_barcode',
                                     'stop_codon',
                                     'too_short',
                                     'total_translated_sequences',
                                     'uag',
                                     'no_left',
                                     'no_right',
                                     'left_and_right',
                                     'not_divisible_by_3',
                                     'divisible_by_3_not_NNK'], 0)
        barcode2statistics[barcode]['lib_type'] = {}

        barcode2info_filenames[barcode] = f'{subdir_path}/{sample_name}_info.txt'

    return barcode2filepaths, barcode2filehandlers, barcode2statistics, barcode2info_filenames


def str_diff(s1, s2):
    # s1 should be equal or shorter. If not, swap.
    if len(s2) < len(s1):
        s1, s2 = s2, s1
    diff = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            diff += 1
    return diff


def write_header(f_handler, txt):
    f_handler.write('_' * 80 + f'\n{txt}')


def filter_reads(argv, fastq_path, parsed_fastq_results, logs_dir,
                 done_path, barcode2samplename_path,
                 left_construct, right_construct, max_mismatches_allowed,
                 min_sequencing_quality, minimal_length_required, gz):

    start_time = datetime.datetime.now()
    from auxiliaries.pipeline_auxiliaries import nnk_table, codon_table

    logger.info(f'{datetime.datetime.now()}: Loading barcode2samplename file from:\n{barcode2samplename_path}')

    if os.path.isdir(fastq_path):
        fastq_files = [os.path.join(fastq_path, file) for file in os.listdir(fastq_path) if file.endswith('.fastq.gz')]
        total_files = len(fastq_files)
        logger.info(f'{total_files} fastq files found: {fastq_files}')
    else:
        fastq_files = [fastq_path]
        total_files = 1

    barcode2samplename = load_table_to_dict(barcode2samplename_path, 'Barcode {} belongs to more than one sample!!')
    assert len(barcode2samplename) > 0, f'No barcodes were found in {barcode2samplename_path}'  # TODO: add informative error to log

    # get barcode's length by taking the length of an arbitrary barcode
    for barcode in barcode2samplename:
        barcode_len = len(barcode)
        break

    # initialize dictionaries that map barcode to output file names, file handlers and statistics
    barcode2filepaths, barcode2filehandlers, barcode2statistics, barcode2info_filenames = get_barcodes_dictionaries(barcode2samplename, parsed_fastq_results, gz)

    lib_types = set()  # keeps all the library types seen
    record_num = 0  # so the line number will represent the actual sequence inside each record
    for i, fastq_file in enumerate(fastq_files):
        logger.info(f'Reading and filtering {i+1}/{total_files}: {fastq_file}')
        with gzip.open(fastq_file, 'rb') as fastq_f:
            for line1 in fastq_f:
                # ignore first line of each record
                line1 = line1.decode("utf-8").rstrip()
                # second line contains the read itself
                dna_read = fastq_f.readline().decode("utf-8").rstrip()
                # ignore third line of each record
                line3 = fastq_f.readline().decode("utf-8").rstrip()
                # fourth line contains the sequencing quality
                quality = fastq_f.readline().decode("utf-8").rstrip()


                if record_num % 1_000_000 == 0:
                    logger.info(f'{datetime.datetime.now()}: {record_num} records have been processed...')

                if not line1.startswith('@'):
                    raise TypeError('NGS file has an invalid format')

                record_num += 1

                # check that barcode exists
                barcode = dna_read[:barcode_len]
                if barcode not in barcode2samplename:
                    continue
                #if barcode is legal
                barcode2statistics[barcode]['legal_barcode'] += 1
                
                # start with the seq and if will be a problem, write it
                barcode2filehandlers[barcode]['filtration_log'].write(f"\nSequence number {barcode2statistics[barcode]['legal_barcode']}: {dna_read}\t") 
                dna = dna_read[8:]
                # verify left construct
                # current_left_construct = dna_read[barcode_len: barcode_len + left_construct_length]
                left_search = regex.search(f'({left_construct}){{s<=1}}', dna, regex.BESTMATCH)
                if left_search == None: # Don't find left
                    barcode2statistics[barcode]['no_left'] += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"left construct has more than 1 mistake\t")
                    continue                   
                else: # Find the right construct
                    right_serach = regex.search(f'({right_construct}){{s<=1}}', dna, regex.BESTMATCH)
                    if right_serach == None:
                        barcode2statistics[barcode]['no_right'] += 1
                        barcode2filehandlers[barcode]['filtration_log'].write(f"right construct has more than 1 mistake\t")
                        continue
                    else: # Left and right!
                        barcode2statistics[barcode]['left_and_right'] += 1
                        


                # extract random fragment (according to the currently assumed library)
                random_dna_start = left_search.span()[1]
                random_dna_end = right_serach.span()[0]
                rest_of_read = dna[random_dna_start:random_dna_end]
                if len(rest_of_read)%3 != 0 :
                    barcode2statistics[barcode]['not_divisible_by_3'] += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"not_divisible_by_3\t")
                    continue
                # for sure divisible by 3
                random_peptide = ''
                has_stop_codon = False
                q_in_peptide = False
                i = 0
                for i in range(0, len(rest_of_read), 3):
                    codon = rest_of_read[i: i+3]
                    if codon == "TGA" or codon == "TAA":
                      has_stop_codon = True
                      barcode2statistics[barcode]['stop_codon'] += 1
                      barcode2filehandlers[barcode]['filtration_log'].write(f"contains a stop codon\t")
                      break

                    if codon[0] not in 'ACGT' or codon[1] not in 'ACGT':  # unrecognized dna base, e.g., N
                        barcode2statistics[barcode]['divisible_by_3_not_NNK'] += 1
                        barcode2filehandlers[barcode]['filtration_log'].write(f"divisible_by_3_not_NNK\t")
                        break
                    if codon[2] not in 'GT':  # K group (of NNK) = G or T
                        barcode2statistics[barcode]['divisible_by_3_not_NNK'] += 1
                        barcode2filehandlers[barcode]['filtration_log'].write(f"divisible_by_3_not_NNK\t")
                        break
                    random_peptide += codon_table[codon]
                    if codon == 'TAG':
                        q_in_peptide = True



                if (len(random_peptide) < minimal_length_required or
                    (random_peptide.startswith('C') and random_peptide.endswith('C') and
                    len(random_peptide)-2 < minimal_length_required)):  # minimum required length (excluding flanking Cysteine)
                    barcode2statistics[barcode]['too_short'] += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"random peptide is too short\t{random_peptide}\t")
                    # all set and documented. We can continue to the next read...
                    continue

                if has_stop_codon:
                    continue

                # reached here? read is valid!! woohoo
                barcode2statistics[barcode]['total_translated_sequences'] += 1
                if q_in_peptide:
                    barcode2statistics[barcode]['uag'] += 1

                # update library counts
                current_lib_type = f'{len(random_peptide)}'
                if random_peptide.startswith('C') and random_peptide.endswith('C'):
                    # A random peptide from a cys loop library
                    current_lib_type = f'C{len(random_peptide)-2}C'
                if current_lib_type not in lib_types:
                    lib_types.add(current_lib_type)
                barcode2statistics[barcode]['lib_type'][current_lib_type] = barcode2statistics[barcode]['lib_type'].get(current_lib_type, 0) + 1


    for barcode in barcode2samplename:    
        barcode2filehandlers[barcode]['filtration_log'].close()

    lib_types = sorted(lib_types, key=lambda x: int(x) if x.isdigit() else 100 + int(x[1:-1]))  # no C's then C's

    for barcode in barcode2samplename:
        with open(barcode2info_filenames[barcode], 'w') as f:
            f.write(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}\n')
            f.write(f'reads_filtration function is invoked with the following parameters:\n')
            f.write(f'fastq_path = {fastq_path}\n'
                    f'parsed_fastq_results = {parsed_fastq_results}\n'
                    f'logs_dir = {logs_dir}\n'
                    f'barcode2samplename_path = {barcode2samplename_path}\n'
                    f'left_construct = {left_construct}\n'
                    f'right_construct = {right_construct}\n'
                    f'max_mismatches_allowed = {max_mismatches_allowed}\n'
                    f'min_sequencing_quality = {min_sequencing_quality}\n')

            total_filtered_reads_per_barcode = barcode2statistics[barcode]['stop_codon'] + \
                                   barcode2statistics[barcode]['too_short'] + \
                                   barcode2statistics[barcode]['no_left'] + \
                                   barcode2statistics[barcode]['no_right'] + \
                                   barcode2statistics[barcode]['not_divisible_by_3'] + \
                                   barcode2statistics[barcode]['divisible_by_3_not_NNK']

            # sanity_check = barcode2statistics[barcode]['poor_quality_barcode'] + barcode2statistics[barcode]['high_quality_barcode'] == barcode2statistics[barcode]['legal_barcode']
            # assert sanity_check, 'poor_quality + high_quality != total'

            write_header(f, f"Total sequences with a appropriate barcode: {barcode2statistics[barcode]['legal_barcode']}\n")
            write_header(f, f"Total number of reads that were filtered (due to the following reasons) -> {total_filtered_reads_per_barcode}\n")
            f.write(f"More than {max_mismatches_allowed} mistakes in the left constant sequences -> {barcode2statistics[barcode]['no_left']}\n"
                    f"More than {max_mismatches_allowed} mistakes in the right constant sequences -> {barcode2statistics[barcode]['no_right']}\n"
                    f"Sequences with left and right constant -> {barcode2statistics[barcode]['left_and_right']}\n"
                    f"Nonsense stop codon -> {barcode2statistics[barcode]['stop_codon']}\n"           
                    f"Too short NNK sequences -> {barcode2statistics[barcode]['too_short']}\n"
                    f"Not divisible by 3 -> {barcode2statistics[barcode]['not_divisible_by_3']}\n"
                    f"Divisible by 3 but not NNK -> {barcode2statistics[barcode]['divisible_by_3_not_NNK']}\n")

            write_header(f, f"\n\nTOTAL NUMBER OF TRANSLATED SEQUENCES -> {barcode2statistics[barcode]['total_translated_sequences']}\n\n")

            write_header(f, 'Translated sequences (per library):\n')
            for lib_type in lib_types:
                f.write(f"{lib_type} -> {barcode2statistics[barcode]['lib_type'].get(lib_type, 0)}\n")

            write_header(f, f"Total number of dna sequences with UAG (q) -> {barcode2statistics[barcode]['uag']}\n")


    logger.info('Writing execution summary file...')

    with open(f'{parsed_fastq_results}/summary_log.txt', 'w') as log_f:
        log_f.write(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}\n')
        log_f.write(f'reads_filtration function is invoked with the following parameters:\n')
        log_f.write(f'fastq_path = {fastq_path}\n'
                    f'parsed_fastq_results = {parsed_fastq_results}\n'
                    f'logs_dir = {logs_dir}\n'
                    f'barcode2samplename_path = {barcode2samplename_path}\n'
                    f'left_construct = {left_construct}\n'
                    f'right_construct = {right_construct}\n'
                    f'max_mismatches_allowed = {max_mismatches_allowed}\n'
                    f'min_sequencing_quality = {min_sequencing_quality}\n')
        log_f.write(f'Barcodes = {" ".join(barcode2samplename)}\n'
                    f'Sample names = {" ".join(barcode2samplename[barcode] for barcode in barcode2samplename)}\n\n')

        reads_with_legal_barcode = sum(barcode2statistics[barcode]['legal_barcode'] for barcode in barcode2samplename)
        write_header(log_f, f'Total number of dna sequences (i.e., fastq records) processed -> {record_num}\n')
        write_header(log_f, f'Total number of dna sequences with an INAPPROPRIATE barcode -> {record_num - reads_with_legal_barcode}\n')
        write_header(log_f, f'Total number of dna sequences with a LEGAL barcode -> {reads_with_legal_barcode}\n')

        for barcode in barcode2samplename:
            log_f.write(f'{barcode} -> {barcode2statistics[barcode]["legal_barcode"]}\n')

        total_filtered_reads = sum(barcode2statistics[barcode]['no_left'] + 
                                barcode2statistics[barcode]['no_right'] + 
                                barcode2statistics[barcode]['not_divisible_by_3'] + 
                                barcode2statistics[barcode]['divisible_by_3_not_NNK'] +
                                barcode2statistics[barcode]['stop_codon'] +
                                barcode2statistics[barcode]['too_short'] for barcode in barcode2samplename)

        write_header(log_f, f"Total number of reads that were filtered (due to the following reasons) -> {total_filtered_reads}\n")

        log_f.write(f"More than {max_mismatches_allowed} mistakes in the left constant sequences -> {sum(barcode2statistics[barcode]['no_left'] for barcode in barcode2samplename)}\n"
                    f"More than {max_mismatches_allowed} mistakes in the right constant sequences -> {sum(barcode2statistics[barcode]['no_right'] for barcode in barcode2samplename)}\n"
                    f"Nonsense stop codon -> {sum(barcode2statistics[barcode]['stop_codon'] for barcode in barcode2samplename)}\n"          
                    f"Too short NNK sequences -> {sum(barcode2statistics[barcode]['too_short'] for barcode in barcode2samplename)}\n"
                    f"Not divisible by 3 -> {sum(barcode2statistics[barcode]['not_divisible_by_3'] for barcode in barcode2samplename)}\n"
                    f"Divisible by 3 not NNK -> {sum(barcode2statistics[barcode]['divisible_by_3_not_NNK'] for barcode in barcode2samplename)}\n")
        
                
        total_translated_sequences = sum(barcode2statistics[barcode]['total_translated_sequences'] for barcode in barcode2samplename)

        write_header(log_f, f'\n\nTOTAL NUMBER OF TRANSLATED SEQUENCES -> {total_translated_sequences}\n\n')
        write_header(log_f, 'No left constant (per barcode):\n')
        for barcode in barcode2samplename:
            log_f.write(f'{barcode} -> {barcode2statistics[barcode]["no_left"]}\n')
        
        write_header(log_f, 'No right constant (per barcode):\n')
        for barcode in barcode2samplename:
            log_f.write(f'{barcode} -> {barcode2statistics[barcode]["no_right"]}\n')
        
        write_header(log_f, 'Translated sequences (per library):\n')
        for lib_type in lib_types:
            log_f.write(f"{lib_type} -> {sum(barcode2statistics[barcode]['lib_type'].get(lib_type, 0) for barcode in barcode2samplename)}\n")

        write_header(log_f, 'Translated peptides (per barcode):\n')
        for barcode in barcode2samplename:
            log_f.write(f'{barcode} -> {barcode2statistics[barcode]["total_translated_sequences"]}\n')

        uag = sum(barcode2statistics[barcode]["uag"] for barcode in barcode2samplename)
        write_header(log_f, f'Total number of dna sequences with UAG (q) -> {uag}\n')

        end_time = datetime.datetime.now()
        write_header(log_f, f'Total running time: {str(end_time-start_time)[:-3]}')

        log_f.close()

        logger.info(f'Started processing at {start_time}')
        logger.info(f'Done Processing! at {end_time}')
        logger.info(f'Total running time: {str(end_time-start_time)[:-3]}')

    with open(done_path, 'w') as f:
        f.write(" ".join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_path', type=str, help='A fastq file to parse')
    parser.add_argument('parsed_fastq_results', type=str, help='folder output')
    parser.add_argument('logs_dir', type=str, help='logs folder')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('barcode2samplename', type=str, help='A path to the barcode to sample name file')
    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('--left_construct', type=str, default="CAACGTGGC", help='left constant sequence')
    parser.add_argument('--right_construct', type=str, default="GCCT", help='right constant sequence')
    parser.add_argument('--max_mismatches_allowed', type=int, default=1,
                        help='number of mismatches allowed together in both constructs')
    parser.add_argument('--min_sequencing_quality', type=int, default=38,
                        help='Minimum average sequencing threshold allowed after filtration'
                             'for more details, see: https://en.wikipedia.org/wiki/Phred_quality_score')
    # more than 12 aa-long random peptide (e.g., C12C) is irrelevant for 51bps-long reads
    # parser.add_argument('--lib_types', type=str.upper, default='6,C6C,8,C8C,10,C10C,12', help='OBSOLETE: Ignore this param. CxC,x')
    parser.add_argument('--minimal_length_required', default=3, type=int,
                        help='Shorter peptides will be discarded')
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

    # try:
    filter_reads(sys.argv, args.fastq_path, args.parsed_fastq_results, args.logs_dir,
                 args.done_file_path, args.barcode2samplename,
                 args.left_construct, args.right_construct, args.max_mismatches_allowed,
                 args.min_sequencing_quality, args.minimal_length_required,
                 True if args.gz else False)
    # except Exception as e:
    #     fail(error_path, e)

