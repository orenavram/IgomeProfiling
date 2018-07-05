#in order to submit this program as a job to the queue, create a copy of the file:
# /groups/pupko/orenavr2/gershoni/Experiments/Exp30/FilterReads.sge and edit it (change the experiment number) accordingly

# The structure is:
# ===============================================
# Barcodes (8 bp) - no mistakes allowed
# fth1 anealed Site (16 bp) - with mistakes allowed
# Domain BC (12 bp)
# fth1 anealed - anti sense (17 bp) - no mistakes allowed


from Auxilaries import *
start = time()

#from subprocess import check_output, call

logger = logging.getLogger('main')
logging.basicConfig(level=logging.INFO) # set level to logging.DEBUG to see debugging comments

if len(argv) < 6:
    logger.error('Usage: python '+argv[0]+' <fastq_path> <out_dir> <barcode_to_name (tab delimited)> <code_to_barcode (tab delimited)> <mistakes_allowed>')
    exit()
    #fastq_path = '/Users/Oren/Dropbox/Projects/gershoni/in_test/first_phase/Exp30.For_DS_Tests_barcodes_GGATC_GGCTA_GGTCA.fastq'
    #out_dir = '/Users/Oren/Dropbox/Projects/gershoni/out_test'
else:
    fastq_path = argv[1]
    out_dir = argv[2].rstrip('/')
    barcode_to_name_file = argv[3]
    code_to_barcode_file = argv[4]
    mistakes_allowed = argv[5]

if len(argv)>6:
    SUFFIX = argv[6]

logger.info('argv: {} ; len(argv)= {}'.format(argv, len(argv)))

def parse_fastq_file(fastq_file, counters, sequence_without_barcode_to_counts):
    local_num_of_seqs_read = 0
    local_num_of_no_barcode = 0
    if fastq_file.endswith('gz'):
        open_method = gzip.open
    else:
        open_method = open
    with open_method(fastq_file) as fastq:
        while True:
            # every read is represented by a text chunk of 4 lines.
            # The sequence is in the second sequence, e.g.:
            # '@D00257:291:CB65BANXX:4:2310:3171:2243 1:N:0:1\nGGATCAAGTAGGGGATCCAGGTGGTGCGCGACCTCTAGAGCCGACCGCGATAGATCGGAAG\n+\nCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n'

            # read next 4-lines chunk
            line1 = fastq.readline().rstrip()
            line2 = fastq.readline().rstrip()
            line3 = fastq.readline().rstrip()
            line4 = fastq.readline().rstrip()

            if not line1:
                break  # end of file

            sequence = line2.rstrip()
            barcode = sequence[:SAMPLE_BARCODE_LENGTH]
            local_num_of_seqs_read += 1

            if local_num_of_seqs_read % 10000 == 0:
                logger.info('num_of_seqs_read is {}'.format(local_num_of_seqs_read))

            if barcode not in barcode_to_name:
                local_num_of_no_barcode += 1
                sequence_without_barcode_to_counts[sequence] = sequence_without_barcode_to_counts.get(sequence, 0) + 1
                continue  # no barcode. skipping sequence...

            # sample barcode exists.
            counters[barcode][TOTAL_SEQ] += 1
            putative_fth1_anealed_site = sequence[SAMPLE_BARCODE_LENGTH:SAMPLE_BARCODE_LENGTH + FTH1_ANEALED_SITE_LENGTH]
            if not regex.match(fth1_anealed_site_regex, putative_fth1_anealed_site):
                counters[barcode][TOO_MANY_MISTAKES] += 1
                N_count = putative_fth1_anealed_site.count('N')  # reflects the sequencing quality
                if N_count:  # unrecognized nucleotides
                    counters[barcode][WRONG_SEQ_WITH_N] += 1
                files_handlers[barcode][FILTERED].write('\t'.join([TOO_MANY_MISTAKES, sequence, 'N:', str(N_count), 'START MISTAKES',
                     str(diff(FTH1_ANEALED_SITE, putative_fth1_anealed_site))]) + '\n')
                continue  # too many mismatches

            # the anealed_site (after the barcode) is as expected or with allowed mistakes
            putative_fth1_anealed_antisense = sequence[ANTISENSE_START:ANTISENSE_START + FTH1_ANEALED_ANTISENSE_LENGTH]
            if not FTH1_ANEALED_ANTISENSE == putative_fth1_anealed_antisense:
                counters[barcode][SEQ_WITH_ERRORS_IN_FTH1_ANNEALED_ANTI_SENSE] += 1
                N_count = putative_fth1_anealed_antisense.count('N')  # reflects the sequencing quality
                if N_count:
                    counters[barcode][WRONG_SEQ_WITH_N] += 1
                files_handlers[barcode][FILTERED].write('\t'.join(['NO_FTH1_Anealed_AntiSense', sequence, 'N:', str(N_count)]) + '\n')
                continue  # Not a PERFECT match

            # the anealed_antisense matches exactly
            # all tests are OK
            counters[barcode][TOTAL_SEQS_OK] += 1
            DomainSeq = sequence[VARIBLE_REGION_START:ANTISENSE_START]
            # N_count = DomainSeq.count('N')  # reflects the sequencing quality
            # files_handlers[barcode][FS].write('_'.join(['>Seq', str(counters[barcode][TOTAL_SEQ]), 'LIB_DS']) + '\n')
            files_handlers[barcode][FS].write(DomainSeq + '\n')
            files_handlers[barcode][FASTQ].write('\n'.join([line1, line2, line3, line4]) + '\n')
    return local_num_of_seqs_read, local_num_of_no_barcode

done_file = os.path.join(out_dir,'done.txt')
if os.path.exists(done_file):
    logger.info('Results for Exp already exist at:\n' + out_dir)
    exit()

#default params
#mistakes_allowed = '1' if len(argv)<4 else argv[3]
#barcode_to_name_file = '/Users/Oren/Dropbox/Projects/gershoni/in_test/first_phase/DS_Samples.txt' if len(argv)<5 else argv[4]
#code_to_barcode_file = '/Users/Oren/Dropbox/Projects/gershoni/in_test/first_phase/code_to_barcode.txt' if len(argv)<6 else argv[5]
#barcode_to_name_file = '/groups/pupko/orenavr2/gershoni/data/EXP_30_DS_Samples.txt' if len(argv)<5 else argv[4]
#code_to_barcode_file = '/groups/pupko/orenavr2/gershoni/data/code_to_barcode.txt' if len(argv)<6 else argv[5]


counter_types = [TOO_MANY_MISTAKES, SEQ_WITH_ERRORS_IN_FTH1_ANNEALED_ANTI_SENSE, WRONG_SEQ_WITH_N, TOTAL_SEQ, TOTAL_SEQS_OK]
#OUT_FILES_TYPES = dict((x, '.'.join(['MistakesAllowed', mistakes_allowed, x])) for x in [FILTERED, FASTQ, INFO, FS])

fth1_anealed_site_regex = '(' + FTH1_ANEALED_SITE + '){s<=' + mistakes_allowed + '}'

barcode_to_name = load_table_to_dict(barcode_to_name_file)
code_to_barcode = load_table_to_dict(code_to_barcode_file)
barcode_to_code = dict([item[::-1] for item in code_to_barcode.items()])

# samples barcodes in the same order at the DS_samples file
barcodes = get_column_from_file(barcode_to_name_file, 0)

logger.debug(barcode_to_code)

files_paths = {}
files_handlers = {}

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

num_of_opened_files = 0

logger.info('num of barcodes is {}. 5 files will be opened for each barcode (4 now + 1 at the end).'.format(len(barcodes)))
#initializing resources
for barcode in barcodes:

    code = barcode_to_code.get(barcode, barcode)
    name = barcode_to_name.get(barcode, barcode)
    logger.info('Barcode: {}'.format(barcode))

    if code == barcode:
        logger.warning('No code for barcode')

    # create sub directories for each sample barcode
    sub_dir = os.path.join(out_dir, name)
    if not os.path.exists(sub_dir):
        os.makedirs(sub_dir)

    # initializing file names and handlers for each out_file_type
    files_paths[barcode] = {}
    files_handlers[barcode] = {}

    #file for domain barcodes that passed the QA
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, FS])) + SUFFIX
    files_paths[barcode][FS] = out_file
    files_handlers[barcode][FS] = open(out_file, 'w')

    #file for filtered sequences
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, FILTERED])) + SUFFIX
    files_paths[barcode][FILTERED] = out_file
    files_handlers[barcode][FILTERED] = open(out_file, 'w')

    #file for filtered + domain barcodes
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, FASTQ])) + '.gz'
    files_paths[barcode][FASTQ] = out_file
    files_handlers[barcode][FASTQ] = gzip.open(out_file, 'w') # these are very (!) large files. Better be kept zipped..

    #info file
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, INFO])) + SUFFIX
    files_paths[barcode][INFO] = out_file
    files_handlers[barcode][INFO] = open(out_file, 'w')

    #counts file
    #NO need to open it now.
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, COUNTS])) + SUFFIX
    files_paths[barcode][COUNTS] = out_file

    num_of_opened_files += 4
    logger.info('{} were opened.'.format(num_of_opened_files))

    if len(files_handlers) != num_of_opened_files/4:
        logger.error('\n' + '#'*50 + '\nThis barcode was already seen!!! Sample barcodes are not unique!!\n' + '#'*50)
        exit()

# counters
num_of_no_barcode = 0
num_of_seqs_read = 0
counters = {}
summary_counters = {}

# per barcode counters
for barcode in barcodes:
    counters[barcode] = {}
    for counter in counter_types:
        counters[barcode][counter] = 0

for counter in counter_types:
    summary_counters[counter] = 0

logger.info('Done initializing parameters! Starting to parse the fastq file...')


#parsing the deep sequencing file
fastq_files = []
if os.path.isdir(fastq_path):
    fastq_files.extend(os.path.join(fastq_path, file_name) for file_name in os.listdir(fastq_path) if 'fastq' in file_name)
else: #fastq_path is a filescp -v
    fastq_files.append(fastq_path)

sequence_without_barcode_to_counts = {}
for fastq_file in fastq_files:
    logger.info('\n' + '#'*50 + '\nStart parsing {}\n'.format(fastq_file) + '#'*50)
    local_values = parse_fastq_file(fastq_file, counters, sequence_without_barcode_to_counts)
    num_of_seqs_read += local_values[0]
    num_of_no_barcode += local_values[1]

logger.info('Done parsing fastq file!')

sorted_sequences_by_counts = sorted(sequence_without_barcode_to_counts, key=sequence_without_barcode_to_counts.get, reverse=True)
text_to_no_barcodes_file = ''
for sequence in sorted_sequences_by_counts:
    text_to_no_barcodes_file += sequence + '\t' + str(sequence_without_barcode_to_counts[sequence]) + '\n'

with open(os.path.join(out_dir, 'no_barcode.txt'), 'w') as f:
    f.write(text_to_no_barcodes_file)

for barcode in barcodes:

    text_to_info_file = ' '.join([argv[0], fastq_path, out_dir, mistakes_allowed]) + '\n'
    text_to_info_file += 'IN: ' + fastq_path + '\n'
    text_to_info_file += 'OUT: '+out_dir+'\n'
    text_to_info_file += 'MISTAKES ALLOWED: '+mistakes_allowed+'\n'
    text_to_info_file += 'BARCODE: '+barcode+'\n'
    text_to_info_file += 'Total seq: '+str(counters[barcode][TOTAL_SEQ]) + '\n'
    text_to_info_file += '\n\nWRONG SEQ STATISTICS\n=============================================\n'
    text_to_info_file += 'Sequences with more than '+mistakes_allowed+' mistakes in fth1 annealed site sequence region: ' + str(counters[barcode][TOO_MANY_MISTAKES]) + '\n'
    text_to_info_file += 'Wrong seq with N: ' + str(counters[barcode][WRONG_SEQ_WITH_N]) + '\n'
    text_to_info_file += 'Seq with errors in the fth1 annealed anti sense construct: ' + str(counters[barcode][SEQ_WITH_ERRORS_IN_FTH1_ANNEALED_ANTI_SENSE]) + '\n'
    text_to_info_file += TOTAL_SEQ_AFTER_FILTERS + ': ' + str(counters[barcode][TOTAL_SEQS_OK]) + '\n'
    files_handlers[barcode][INFO].write(text_to_info_file)

    for counter in counter_types:
        summary_counters[counter] += counters[barcode][counter]


text_to_info_file = ' '.join([argv[0], fastq_path, out_dir, mistakes_allowed]) + '\n'
text_to_info_file += 'IN: ' + fastq_path + '\n'
text_to_info_file += 'OUT: '+out_dir+'\n'
text_to_info_file += 'MISTAKES ALLOWED: '+mistakes_allowed+'\n'
text_to_info_file += 'BARCODES: '+','.join(barcode_to_name)+'\n'
text_to_info_file += '\n\n=============================================\n'
text_to_info_file += 'TOTAL ANALYZED SEQ: ' + str(num_of_seqs_read) + '\n'
text_to_info_file += 'SKIPPED SEQ WITH NO BARCODE: ' + str(num_of_no_barcode) + '\n'

for barcode in barcodes:
    text_to_info_file += barcode + ': ' + str(counters[barcode][TOTAL_SEQ]) + '\n'

text_to_info_file += '\nTotal seq with barcodes: ' + str(summary_counters[TOTAL_SEQ]) + '\n'
text_to_info_file += '\n\nWRONG SEQ STATISTICS\n=============================================\n'
text_to_info_file += 'Sequences with more than '+mistakes_allowed+' mistakes in fth1 annealed site sequence region: ' + str(summary_counters[TOO_MANY_MISTAKES]) + '\n'
text_to_info_file += 'Wrong seq with N: ' + str(summary_counters[WRONG_SEQ_WITH_N]) + '\n'
text_to_info_file += 'Seq with no stop codon construct: ' + str(summary_counters[SEQ_WITH_ERRORS_IN_FTH1_ANNEALED_ANTI_SENSE]) + '\n'
text_to_info_file += '\n\nTotal sequences after filters - by barcode:\n=============================================\n'

for barcode in barcodes:
    text_to_info_file += barcode + ': ' + str(counters[barcode][TOTAL_SEQS_OK]) + '\n'

text_to_info_file += TOTAL_SEQ_AFTER_FILTERS + ': ' + str(summary_counters[TOTAL_SEQS_OK]) + '\n'

with open(os.path.join(out_dir, 'info.txt'), 'w') as f:
    f.write(text_to_info_file)

logger.info('Freeing up resources...')
for barcode in files_handlers:
    for handler in files_handlers[barcode]:
        files_handlers[barcode][handler].close()

#frequency_counters_per_sample = {}
logger.info('Writing unique files...')
for barcode in barcodes:
    logger.info('Barcode: ' + barcode)
    #frequency_counters_per_sample[barcode] = {}
    with open(files_paths[barcode][FS]) as f:
        txt = f.read()
    with open(files_paths[barcode][COUNTS], 'w') as f:
        counts = ''
        frequency_counter = Counter(txt.split())
        # with open(files_paths[barcode][COUNTS][:-4]+'2'+files_paths[barcode][COUNTS][-4:]) as f2:
        #     txt = f2.read()
        #     lines = [x.split('_')[-1].strip().split()[::-1] for x in txt.split('>')[1:]]
        #     d = Counter(dict([(x, int(y)) for x,y in lines]))
        #     print d.most_common()
        #     print frequency_counter.most_common()
        #     assert d==frequency_counter
        for seq, count in frequency_counter.most_common():
            counts += '\t'.join([seq, str(count)]) + '\n'
            #frequency_counters_per_sample[barcode][seq] = count
        f.write(counts)
    os.remove(files_paths[barcode][FS])

#touch done
with open(done_file, 'w') as f:
    pass

end = time()
logger.info('Finished filtering reads. Took {}'.format(measure_time(start, end)))


