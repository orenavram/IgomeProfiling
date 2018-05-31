#in order to submit this program as a job to the queue, create a copy of the file:
# /groups/pupko/orenavr2/gershoni/Experiments/Exp30/FilterReads.sge and edit it (change the experiment number) accordingly

# The structure is:
# ===============================================
# Barcodes (8 bp) - no mistakes allowed
# fth1 anealed Site (16 bp) - with mistakes allowed
# Domain BC (12 bp)
# fth1 anealed - anti sense (17 bp) - no mistakes allowed

from time import time
start = time()

#from subprocess import check_output, call
from collections import Counter
from Auxilaries import *
import logging
logger = logging.getLogger('main')
logging.basicConfig(level=logging.INFO) # set level to logging.DEBUG to see debugging comments

if len(argv) < 6:
    logger.error('Usage: python '+argv[0]+' <fastq_file> <out_dir> <barcode_to_name (tab delimited)> <code_to_barcode (tab delimited)> <mistakes_allowed>')
    exit()
    #fastq_file = '/Users/Oren/Dropbox/Projects/gershoni/in_test/first_phase/Exp30.For_DS_Tests_barcodes_GGATC_GGCTA_GGTCA.fastq'
    #out_dir = '/Users/Oren/Dropbox/Projects/gershoni/out_test'
else:
    fastq_file = argv[1]
    out_dir = argv[2].rstrip('/')
    barcode_to_name_file = argv[3]
    code_to_barcode_file = argv[4]
    mistakes_allowed = argv[5]

if len(argv)>6:
    suffix = argv[6]

logger.info('argv: {} ; len(argv)= {}'.format(argv, len(argv)))

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


counter_types = [TooManyMistakes, Seq_with_errors_in_fth1_annealed_anti_sense, Wrong_Seq_With_N, TotalSeq, TotalSeqsOK]
#OUT_FILES_TYPES = dict((x, '.'.join(['MistakesAllowed', mistakes_allowed, x])) for x in [FILTERED, FASTQ, INFO, FS])

fth1_anealed_site_regex = '('+fth1_anealed_site+'){s<='+mistakes_allowed+'}'

barcode_to_name = load_table_to_dict(barcode_to_name_file)
code_to_barcode = load_table_to_dict(code_to_barcode_file)
barcode_to_code = dict([item[::-1] for item in code_to_barcode.items()])

logger.debug(barcode_to_code)

files_paths = {}
files_handlers = {}

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

#initializing resources
for barcode in barcode_to_name:

    code = barcode_to_code.get(barcode, barcode)
    name = barcode_to_name.get(barcode, barcode)
    if code == barcode:
        logger.warning('\n'+('!'*5)+'Warning: No code/name for barcode'+barcode+('!'*5))

    # create sub directories for each sample barcode
    sub_dir = os.path.join(out_dir, name)
    if not os.path.exists(sub_dir):
        os.makedirs(sub_dir)

    # initializing file names and handlers for each out_file_type
    files_paths[barcode] = {}
    files_handlers[barcode] = {}

    #file for domain barcodes that passed the QA
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, FS])) + suffix
    files_paths[barcode][FS] = out_file
    files_handlers[barcode][FS] = open(out_file, 'w')

    #file for filtered sequences
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, FILTERED])) + suffix
    files_paths[barcode][FILTERED] = out_file
    files_handlers[barcode][FILTERED] = open(out_file, 'w')

    #file for filtered + domain barcodes
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, FASTQ])) + '.gz'
    files_paths[barcode][FASTQ] = out_file
    files_handlers[barcode][FASTQ] = gzip.open(out_file, 'w') # these are very (!) large files. Better be kept zipped..

    #info file
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, INFO])) + suffix
    files_paths[barcode][INFO] = out_file
    files_handlers[barcode][INFO] = open(out_file, 'w')

    #counts file
    #NO need to open it now.
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', mistakes_allowed, COUNTS])) + suffix
    files_paths[barcode][COUNTS] = out_file
    files_handlers[barcode][COUNTS] = open(out_file, 'w')

# counters
no_barcode = 0
num_of_seqs_read = 0
counters = {}
summary_counters = {}

# per barcode counters
for barcode in barcode_to_name:
    counters[barcode] = {}
    for counter in counter_types:
        counters[barcode][counter] = 0

for counter in counter_types:
    summary_counters[counter] = 0

logger.info('Done initializing parameters! Starting to parse the fastq file...')

#parsing the deep sequencing file
with gzip.open(fastq_file) as fastq:
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
            break #end of file

        sequence = line2.rstrip()
        barcode = sequence[:sample_barcode_length]
        num_of_seqs_read += 1

        if num_of_seqs_read % 10000 == 0:
            logger.info('num_of_seqs_read is {}'.format(num_of_seqs_read))

        if barcode not in barcode_to_name:
            no_barcode += 1
            continue #no barcode. skipping sequnce...

        # sample barcode exists.
        counters[barcode][TotalSeq] += 1
        putative_fth1_anealed_site = sequence[sample_barcode_length:sample_barcode_length+fth1_anealed_site_length]
        if not regex.match(fth1_anealed_site_regex, putative_fth1_anealed_site):
            counters[barcode][TooManyMistakes] += 1
            N_count = putative_fth1_anealed_site.count('N')  # reflects the sequencing quality
            if N_count: #unrecognized nucleotides
                counters[barcode][Wrong_Seq_With_N] += 1
            files_handlers[barcode][FILTERED].write('\t'.join([TooManyMistakes, sequence, 'N:', str(N_count), 'START MISTAKES', str(diff(fth1_anealed_site, putative_fth1_anealed_site))]) + '\n')
            continue #too many mismatches

        # the anealed_site (after the barcode) is as expected or with allowed mistakes
        putative_fth1_anealed_antisense = sequence[antisense_starts:antisense_starts+fth1_anealed_antisense_length]
        if not fth1_anealed_antisense == putative_fth1_anealed_antisense:
            counters[barcode][Seq_with_errors_in_fth1_annealed_anti_sense] += 1
            N_count = putative_fth1_anealed_antisense.count('N')  # reflects the sequencing quality
            if N_count:
                counters[barcode][Wrong_Seq_With_N] += 1
            files_handlers[barcode][FILTERED].write('\t'.join(['NO_FTH1_Anealed_AntiSense', sequence, 'N:', str(N_count)]) + '\n')
            continue #Not a PERFECT match

        # the anealed_antisense matches exactly
        # all tests are OK
        counters[barcode][TotalSeqsOK] += 1
        DomainSeq = sequence[var_starts:antisense_starts]
        #N_count = DomainSeq.count('N')  # reflects the sequencing quality
        #files_handlers[barcode][FS].write('_'.join(['>Seq', str(counters[barcode][TotalSeq]), 'LIB_DS']) + '\n')
        files_handlers[barcode][FS].write(DomainSeq + '\n')
        files_handlers[barcode][FASTQ].write('\n'.join([line1, line2, line3, line4]) + '\n')

logger.info('Done parsing fastq file!')

for barcode in barcode_to_name:

    text_to_info_file = ' '.join([argv[0], fastq_file, out_dir, mistakes_allowed])+'\n'
    text_to_info_file += 'IN: ' + fastq_file + '\n'
    text_to_info_file += 'OUT: '+out_dir+'\n'
    text_to_info_file += 'MISTAKES ALLOWED: '+mistakes_allowed+'\n'
    text_to_info_file += 'BARCODE: '+barcode+'\n'
    text_to_info_file += 'Total seq: '+str(counters[barcode][TotalSeq])+'\n'
    text_to_info_file += '\n\nWRONG SEQ STATISTICS\n=============================================\n'
    text_to_info_file += 'Sequences with more than '+mistakes_allowed+' mistakes in fth1 annealed site sequence region: ' + str(counters[barcode][TooManyMistakes]) + '\n'
    text_to_info_file += 'Wrong seq with N: ' + str(counters[barcode][Wrong_Seq_With_N]) + '\n'
    text_to_info_file += 'Seq with errors in the fth1 annealed anti sense construct: ' + str(counters[barcode][Seq_with_errors_in_fth1_annealed_anti_sense]) + '\n'
    text_to_info_file += TOTAL_SEQ_AFTER_FILTERS + ': ' + str(counters[barcode][TotalSeqsOK]) + '\n'
    files_handlers[barcode][INFO].write(text_to_info_file)

    for counter in counter_types:
        summary_counters[counter] += counters[barcode][counter]


text_to_info_file = ' '.join([argv[0], fastq_file, out_dir, mistakes_allowed])+'\n'
text_to_info_file += 'IN: ' + fastq_file + '\n'
text_to_info_file += 'OUT: '+out_dir+'\n'
text_to_info_file += 'MISTAKES ALLOWED: '+mistakes_allowed+'\n'
text_to_info_file += 'BARCODES: '+','.join(barcode_to_name)+'\n'
text_to_info_file += '\n\n=============================================\n'
text_to_info_file += 'TOTAL ANALYZED SEQ: ' + str(num_of_seqs_read) + '\n'
text_to_info_file += 'SKIPPED SEQ WITH NO BARCODE: ' + str(no_barcode) + '\n'

for barcode in barcode_to_name:
    text_to_info_file += barcode + ': ' + str(counters[barcode][TotalSeq]) + '\n'

text_to_info_file += '\nTotal seq with barcodes: ' + str(summary_counters[TotalSeq]) + '\n'
text_to_info_file += '\n\nWRONG SEQ STATISTICS\n=============================================\n'
text_to_info_file += 'Sequences with more than '+mistakes_allowed+' mistakes in fth1 annealed site sequence region: ' + str(summary_counters[TooManyMistakes]) + '\n'
text_to_info_file += 'Wrong seq with N: ' + str(summary_counters[Wrong_Seq_With_N]) + '\n'
text_to_info_file += 'Seq with no stop codon construct: ' + str(summary_counters[Seq_with_errors_in_fth1_annealed_anti_sense]) + '\n'
text_to_info_file += '\n\nTotal sequences after filters - by barcode:\n=============================================\n'

for barcode in barcode_to_name:
    text_to_info_file += barcode + ': ' + str(counters[barcode][TotalSeqsOK]) + '\n'

text_to_info_file += TOTAL_SEQ_AFTER_FILTERS + ': ' + str(summary_counters[TotalSeqsOK]) + '\n'

with open(os.path.join(out_dir, 'info.txt'), 'w') as f:
    f.write(text_to_info_file)

# freeing up resources
for barcode in files_handlers:
    for handler in files_handlers[barcode]:
        files_handlers[barcode][handler].close()

#frequency_counters_per_sample = {}
#write unique files
for barcode in files_paths:
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

#touch done
with open(done_file, 'w') as f:
    pass

end = time()
logger.info('Finished filtering reads. Took {}'.format(measure_time(start, end)))


