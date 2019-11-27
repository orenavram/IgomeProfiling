import gzip
import regex
import os
import datetime


def load_barcode_to_sample_name(barcode2samplename_path):
    barcode_to_samplename = {}
    with open(barcode2samplename_path) as f:
        for line in f:
            if line.isspace():  # empty line
                continue
            barcode, sample_name = line.strip().split()
            if barcode in barcode_to_samplename:
                assert False, f'Barcode {barcode} belongs to more than one sample!!'  # TODO: write to a global error log file
            barcode_to_samplename[barcode] = sample_name
    return barcode_to_samplename


def get_barcodes_dictionaries(barcode_to_samplename, lib_types, output_dir) -> {str: {str: str}}:
    barcode2filenames = {}
    barcode2filehandlers = {}
    barcode2statistics = {}
    barcode2info_filenames = {}

    for barcode in barcode_to_samplename:
        sample_name = barcode_to_samplename[barcode]  # f'{barcode}_{os.path.split(fastq_file)[-1]}'

        # create a sub dir for that barcode (under output_dir)
        subdir_path = f'{output_dir}/{sample_name}_{barcode}'
        os.makedirs(subdir_path, exist_ok=True)

        # add filenames
        barcode2filenames[barcode] = {}
        # translated seq # TODO: change to .gz
        barcode2filenames[barcode]['fastq'] = f'{subdir_path}/{sample_name}.fastq'
        # translated seq # TODO: change to .gz
        barcode2filenames[barcode]['fna'] = f'{subdir_path}/{sample_name}.fna'
        # translated seq
        barcode2filenames[barcode]['faa'] = f'{subdir_path}/{sample_name}.faa'
        # where in the filter seqs dropped
        barcode2filenames[barcode]['filtration_log'] = f'{subdir_path}/{sample_name}.filtration_log.txt'

        # add filehandlers
        barcode2filehandlers[barcode] = {}
        # TODO: change to zip.open, 'wb'
        barcode2filehandlers[barcode]['fastq'] = open(barcode2filenames[barcode]['fastq'], 'w')
        # TODO: change to zip.open, 'wb'
        barcode2filehandlers[barcode]['fna'] = open(barcode2filenames[barcode]['fna'], 'w')
        barcode2filehandlers[barcode]['faa'] = open(barcode2filenames[barcode]['faa'], 'w')
        barcode2filehandlers[barcode]['filtration_log'] = open(barcode2filenames[barcode]['filtration_log'], 'w')

        barcode2statistics[barcode] = dict.fromkeys(['total_sequences',
                                     'poor_quality_barcode',
                                     'high_quality_barcode',
                                     'too_many_mistakes',
                                     'not_nnk',
                                     'stop_codon',
                                     'total_translated_sequences',
                                     'uag',
                                     'N_count',
                                     'no_library_matched',
                                     'seq_with_stop_codon',
                                     *lib_types], 0)

        barcode2info_filenames[barcode] = f'{subdir_path}/{sample_name}_info.txt'

    return barcode2filenames, barcode2filehandlers, barcode2statistics, barcode2info_filenames


def nnk_translate(dna_sequence, nnk_table):
    aa_seq = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        if codon in nnk_table:
            aa = nnk_table[codon]
        else:
            if codon == "TGA" or codon == "TAA":
                return "stop_codon"
            else:
                return "not_nnk"
        aa_seq += aa
    return aa_seq


def is_cys_loop(seq, nnk_table):
    return seq[:3] in nnk_table and \
           seq[-3:] in nnk_table and\
           nnk_table[seq[:3]] == nnk_table[seq[-3:]] == 'C'


def str_diff(s1, s2):
    # s1 should be equal or shorter. If not, swap.
    if len(s2) < len(s1):
        s1, s2 = s2, s1
    diff = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            diff += 1
    return diff

def filter_reads(argv, fastq_file, output_dir, barcode2samplename_path,
                 left_construct, right_construct, max_mismatches_allowed,
                 min_sequencing_quality, lib_types):

    start_time = datetime.datetime.now()

    nnk_table: {str: str} = {"CGT": "R", "CGG": "R", "AGG": "R",
                 "CTT": "L", "CTG": "L", "TTG": "L",
                 "TCT": "S", "TCG": "S", "AGT": "S",
                 "GCT": "A", "GCG": "A",
                 "GGT": "G", "GGG": "G",
                 "CCT": "P", "CCG": "P",
                 "ACT": "T", "ACG": "T",
                 "CAG": "Q",
                 "TAG": "q",
                 "GTT": "V", "GTG": "V",
                 "AAT": "N",
                 "GAT": "D",
                 "TGT": "C",
                 "GAG": "E",
                 "CAT": "H",
                 "ATT": "I",
                 "AAG": "K",
                 "ATG": "M",
                 "TTT": "F",
                 "TGG": "W",
                 "TAT": "Y"}

    left_construct_length = len(left_construct)
    right_construct_length = len(right_construct)

    logger.info(f'{datetime.datetime.now()}: Loading barcode2samplename file from:\n{barcode2samplename_path}')
    barcode_to_samplename = load_barcode_to_sample_name(barcode2samplename_path)
    assert len(barcode_to_samplename) > 0, f'No barcodes were found in {barcode2samplename_path}'  # TODO: add informative error to log

    # get barcode's length by taking the length of an arbitrary barcode
    for barcode in barcode_to_samplename:
        barcode_len = len(barcode)
        break

    # initialize dictionaries that map barcode to output file names, file handlers and statistics
    barcode2filenames, barcode2filehandlers, barcode2statistics, barcode2info_filenames = get_barcodes_dictionaries(barcode_to_samplename, lib_types, output_dir)

    log_f = open(f'{output_dir}/summary_log.txt', 'w')

    # counters:
    no_barcode = 0
    poor_quality = 0
    too_many_mistakes = 0
    not_nnk = 0
    stop_codon = 0

    record_num = 0  # so the line number will represent the actual sequence inside each record
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
            if barcode not in barcode_to_samplename:
                no_barcode += 1
                continue

            barcode2statistics[barcode]['total_sequences'] += 1

            # TODO: remove
            a = barcode_to_samplename[barcode]
            b = barcode2statistics[barcode]['total_sequences']

            # check that the barcode quality is above the required threshold
            barcode_quality = quality[:barcode_len]
            lowest_scored_character = min([ord(c) for c in barcode_quality])
            if lowest_scored_character < min_sequencing_quality:
                barcode2statistics[barcode]['poor_quality_barcode'] += 1
                poor_quality += 1
                barcode2filehandlers[barcode]['filtration_log'].write(f"{barcode2statistics['total_sequences']}\tlow quality barcode ({lowest_scored_character})\t{dna_read}\n")
                continue
            barcode2statistics[barcode]['high_quality_barcode'] += 1

            # verify left construct
            current_left_construct = dna_read[barcode_len: barcode_len + left_construct_length]
            if regex.match(f'({current_left_construct})' + '{s<='f"{max_mismatches_allowed}"'}', left_construct) == None:
                barcode2statistics[barcode]['too_many_mistakes'] += 1
                too_many_mistakes += 1
                barcode2filehandlers[barcode]['filtration_log'].write(f"{barcode2statistics[barcode]['total_sequences']}\tleft construct has {num_of_mismatches} (which is more than the allowed amount: {max_mismatches_allowed} mismatches)\t{dna_read}\n")
                continue

            # if an appropriate library is found this variable will be set to False
            problematic_sequence = True

            # TODO: make one file for short reads analysis and one for long ones

            for current_lib_type in lib_types:
                # Since the reads are too short (we have 50 bps out of which 5 are barcode's and 9 are left construct's)
                # there are only 36 bps left for the random fragment + 4 bps of the right construct. Thus, we cannot
                # infer directly from which library the read comes from. What we do is the following- for each possible
                # library, we assume it is the correct one and see if the read's configuration (library) makes sense. If
                # it does, we stop searching for "the right configuration". The order of the library types is very
                # important as 2 configurations might fit but we will choose (arbitrarily) the first one.
                current_library_has_c_loop = False
                if current_lib_type[0] == current_lib_type[-1] == 'C':
                    current_lib_aa_length = int(current_lib_type[1:-1])  # e.g., 'C10C' -> 10
                    current_library_has_c_loop = True
                else:
                    current_lib_aa_length = int(current_lib_type)  # e.g., '10' -> 10

                current_library_dna_length = 3 * current_lib_aa_length
                if current_library_has_c_loop:
                    current_library_dna_length += 6 # for the flanking C's

                # extract random fragment (according to the currently assumed library)
                random_dna_start = barcode_len + left_construct_length
                random_dna_end = random_dna_start + current_library_dna_length
                random_dna = dna_read[random_dna_start: random_dna_end]

                if current_library_has_c_loop and not is_cys_loop(random_dna, nnk_table):
                    # library does not fit
                    continue

                # the right construct is induced by the current library
                # not sure if it will success (maybe this library won't fit and we will have to try another one)
                current_right_construct = dna_read[random_dna_end: random_dna_end + right_construct_length]
                num_of_mismatches = str_diff(current_right_construct, right_construct)
                if num_of_mismatches > max_mismatches_allowed:
                    barcode2statistics[barcode]['too_many_mistakes'] += 1
                    too_many_mistakes += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"{barcode2statistics[barcode]['total_sequences']}\tright construct has {num_of_mismatches} (which is more than the allowed amount: {max_mismatches_allowed} mismatches)\t{dna_read}\n")
                    continue

                # maybe it's a C12C but the read is too short so we see C11
                # C12C requires a special handling because it is too long for the short reads (50 bps)
                # 5 (barcode) + 9 (right construct) + 36 (C11) -> 50 (total read's length)
                if current_lib_type == 'C12C':
                    # This block should never happen if "C12C" is checked lastly because the read will be assigned as "12"
                    logger.error(f"{datetime.datetime.now()}: >Seq_{barcode2statistics[barcode]['total_sequences']}_Lib_{current_lib_type}\n{dna_read}\n")
                    logger.error(f"{datetime.datetime.now()}: >Seq_{barcode2statistics[barcode]['total_sequences']}_Lib_{current_lib_type}\n{random_peptide}\n")
                    raise ValueError('How come? "C12C" is checked lastly so this read should have been assigned as "12"')
                    pass

                random_peptide = nnk_translate(random_dna, nnk_table)
                if random_peptide == 'not_nnk':
                    barcode2statistics[barcode]['not_nnk'] += 1
                    not_nnk += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"{barcode2statistics[barcode]['total_sequences']}\trandom peptide is not nnk\t{dna_read[random_dna_start:]}\n")
                elif random_peptide == 'stop_codon':
                    barcode2statistics[barcode]['stop_codon'] += 1
                    stop_codon += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"{barcode2statistics[barcode]['total_sequences']}\tcontains a stop codon\t{dna_read[random_dna_start:]}\n")
                else:
                    # reached here? read is valid!! woohoo
                    barcode2statistics[barcode]['total_translated_sequences'] += 1
                    barcode2statistics[barcode][current_lib_type] += 1
                    # add record to a dedicated fastq file
                    barcode2filehandlers[barcode]['fastq'].write(f'{line1}\n{dna_read}\n{line3}\n{quality}\n')
                    # add dna read to a dedicated fna file
                    barcode2filehandlers[barcode]['fna'].write(f">Seq_{barcode2statistics[barcode]['total_sequences']}_Lib_{current_lib_type}\n{random_dna}\n")
                    # add peptide (translated dna) to a dedicated faa file
                    barcode2filehandlers[barcode]['faa'].write(f">Seq_{barcode2statistics[barcode]['total_sequences']}_Lib_{current_lib_type}\n{random_peptide}\n")

                    if 'q' in random_peptide:
                        barcode2statistics[barcode]['uag'] += 1

                break

            else:
                # no library matched (assuming lib_types is not empty...)
                barcode2statistics[barcode]['no_library_matched'] += 1

                N_count = dna_read[: random_dna_end + right_construct_length].count('N')
                if N_count > 0:
                    barcode2statistics[barcode]['N_count'] += 1

                barcode2filehandlers[barcode]['filtration_log'].write(f"{barcode2statistics[barcode]['total_sequences']}\tno library matched (total N's count is: {N_count})\t{dna_read[random_dna_start:]}\n")


    for barcode in barcode_to_samplename:
        barcode2filehandlers[barcode]['fastq'].close()
        barcode2filehandlers[barcode]['faa'].close()
        barcode2filehandlers[barcode]['fna'].close()
        barcode2filehandlers[barcode]['filtration_log'].close()


    for barcode in barcode_to_samplename:
        with open(barcode2info_filenames[barcode], 'w') as f:
            f.write(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}\n')
            f.write(f'filter_reads function is invoked with the following parameters:\n')
            f.write(f'fastq_file = {fastq_file}\n'
                    f'output_dir = {output_dir}\n'
                    f'barcode2samplename_path = {barcode2samplename_path}\n'
                    f'left_construct = {left_construct}\n'
                    f'right_construct = {right_construct}\n'
                    f'max_mismatches_allowed = {max_mismatches_allowed}\n'
                    f'min_sequencing_quality = {min_sequencing_quality}\n'
                    f'lib_types = {lib_types}\n'
                    f'\n\n\n')
            f.write('')

            sanity_check = barcode2statistics[barcode]['poor_quality_barcode'] + barcode2statistics[barcode]['high_quality_barcode'] == barcode2statistics[barcode]['total_sequences']
            assert sanity_check, 'poor_quality + high_quality != total'
            f.write(f"Total sequences with an appropriate barcode: {barcode2statistics[barcode]['total_sequences']}\n")
            f.write(f'Filtration statistics:\n'
                    f"Poor quality barcode -> {barcode2statistics[barcode]['poor_quality_barcode']}"
                    f"(Sanity check: {sanity_check}"
                    f"{barcode2statistics[barcode]['poor_quality_barcode']} (poor) + "
                    f"{barcode2statistics[barcode]['high_quality_barcode']} (high) = "
                    f"{barcode2statistics[barcode]['total_sequences']} (total)\n"
                    f"Too many mismatches in flanking constant sequences -> {barcode2statistics[barcode]['too_many_mistakes']}\n"
                    f"Not NNK sequences -> {barcode2statistics[barcode]['not_nnk']}\n"
                    f"Nonsense stop codon -> {barcode2statistics[barcode]['stop_codon']}\n"
                    f"UAG stop codon (q) -> {barcode2statistics[barcode]['uag']}\n"
                    f"No library matched -> {barcode2statistics[barcode]['no_library_matched']}\n"
                    f"Too many unrecognized nucleotides (N) -> {barcode2statistics[barcode]['N_count']}\n"
                    f"Total sequences after filtration -> {barcode2statistics[barcode]['total_translated_sequences']}\n"
                    f'\n\n\n')



    log_f.close()
    logger.info(f'Started processing at {start_time}')
    logger.info(f'Done Processing! at {datetime.datetime.now()}')




if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_file', type=str, help='fastq file')
    parser.add_argument('output_dir', type=str, help='folder output')
    parser.add_argument('barcode2samplename', type=str, help='A path to the barcode to sample name file')
    # parser.add_argument('--gz', action='store_true', help='Turn on if the fastq is zipped (False by default)')
    parser.add_argument('--left_construct', type=str, default="CAACGTGGC", help='SfilSite')
    parser.add_argument('--right_construct', type=str, default="GCCT", help='RightConstruct')
    parser.add_argument('--MistakeAllowed', type=int, default=1,
                        help='number of mismatches allowed in EACH construct')
    parser.add_argument('--min_sequencing_quality', type=int, default=38,
                        help='Minimum average sequencing threshold allowed after filtration'
                             'for more details, see: https://en.wikipedia.org/wiki/Phred_quality_score')
    parser.add_argument('--lib_types', type=str.upper, default='C6C,C8C,C10C,6,8,10,12,7,C12C', help='CxC,x')  # However, this makes more sense: 'C12C,C10C,12,C8C,10,C6C,8,6'
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    filter_reads(argv, args.fastq_file, args.output_dir, args.barcode2samplename,
                 args.left_construct, args.right_construct, args.MistakeAllowed,
                 args.min_sequencing_quality, args.lib_types.split(','))
