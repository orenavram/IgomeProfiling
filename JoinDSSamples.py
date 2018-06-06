from time import time
start = time()

from Auxilaries import *

import logging
logger = logging.getLogger('main')
logging.basicConfig(level=logging.INFO) # set level to logging.DEBUG to see debugging comments

def split_dict_by_interest(domains_of_interest, domain_to_counts, correction):
    of_interest = dict.fromkeys(domains_of_interest, '0') # set to 0 by default
    of_interest_corrected = dict.fromkeys(domains_of_interest, correction) # set to $correction by default
    leftovers = {}
    for domain in domain_to_counts:
        if domain in domains_of_interest:
            of_interest[domain] = domain_to_counts[domain]
            of_interest_corrected[domain] = str(int(domain_to_counts[domain]) + int(of_interest_corrected[domain]))
            print(of_interest_corrected[domain])
        else:
            leftovers[domain] = domain_to_counts[domain]
    return of_interest, of_interest_corrected, leftovers

if len(argv) < 3:
    # logger.error('Usage: python '+argv[0]+' <names_to_domains_of_interest_file> <out_dir>; <?correction>')
    # exit(-1)
    names_to_domains_of_interest_file = '/Users/Oren/Dropbox/Projects/gershoni/domains_of_interest/HCV-HIV_P8.txt'
    out_dir = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp27B/first_phase_output/'
    correction = '50'
    lib_name = 'p8'.lower()
else:
    names_to_domains_of_interest_file = argv[1]
    out_dir = argv[2].rstrip('/')
    correction = '50' if len(argv)<=3 else argv[3]
    lib_name = os.path.split(names_to_domains_of_interest_file)[-1].split('_')[1].split('.')[0].lower()

logger.info('argv: {} ; len(argv)= {}'.format(argv, len(argv)))
logger.info('correction={}'.format(correction))
logger.info('lib_name={}'.format(lib_name))
copy_number_file = os.path.join(out_dir, lib_name + '.copy_number.csv')
corrected_copy_number_file = os.path.join(out_dir, lib_name + '.copy_number.corrected_by_'+correction+'.csv')
corrected_percents_file = os.path.join(out_dir, lib_name + '.percent.corrected_by_'+correction+'.csv')
real_vs_leftovers_file = os.path.join(out_dir, lib_name + '.real_vs_leftovers.csv')

#counter_types = [TOO_MANY_MISTAKES, SEQ_WITH_ERRORS_IN_FTH1_ANNEALED_ANTI_SENSE, WRONG_SEQ_WITH_N, TOTAL_SEQ, TOTAL_SEQS_OK]
#OUT_FILES_TYPES = dict((x, '.'.join(['MistakesAllowed', mistakes_allowed, x])) for x in [FILTERED, FASTQ, INFO, FS])


domains_of_interest_to_names = load_table_to_dict(names_to_domains_of_interest_file, backwards=True)
names_to_domains_of_interest = load_table_to_dict(names_to_domains_of_interest_file)
#read the list of relevant domains
domains_of_interest_list = domains_of_interest_to_names.keys()

copy_number_result = ','.join(['sample'] + [domains_of_interest_to_names[domain] for domain in domains_of_interest_list] + ['TOTAL']) + '\n'
corrected_copy_number_result = copy_number_result
corrected_percents_result = ','.join(['sample'] + [domains_of_interest_to_names[domain] for domain in domains_of_interest_list]) + '\n'
real_vs_leftovers_result = ','.join(['real_file', 'sum_real', 'percent_real', 'leftover_file', 'sum_leftover', 'percent_leftover', 'total_reads']) + '\n'

data_is_not_empty = False

logger.info('Aggregating the following files:')
#initializing resources
for sub_dir_name in os.listdir(out_dir):
    sub_dir_path = os.path.join(out_dir, sub_dir_name)
    logger.debug('sub_dir_path ' + sub_dir_path)
    if not lib_name in sub_dir_name.lower():
        logger.debug('lib_name {} NOT in sub_dir_name {}'.format(lib_name, sub_dir_name.lower()))
    if os.path.isdir(sub_dir_path) and lib_name in sub_dir_name.lower(): # scan only relevant dirs
        for file_name in os.listdir(sub_dir_path):
            logger.debug('file_name ' + file_name)
            if file_name.endswith(COUNTS + SUFFIX):
                counts_file_path = os.path.join(sub_dir_path, file_name)
                domain_to_counts = load_table_to_dict(counts_file_path, open_operator=open)
                domain_of_interest_to_counts, domain_of_interest_to_corrected_counts, domain_leftovers_to_counts =\
                    split_dict_by_interest(set(domains_of_interest_list), domain_to_counts, correction)

                logger.info('Aggregating {}'.format(counts_file_path))

                names=[]
                with open(names_to_domains_of_interest_file) as f:
                    for line in f:
                        name = line.split()[0]
                        names.append(name)
                sorted_domain_of_interest_by_name = [names_to_domains_of_interest[name] for name in names]
                print(sorted_domain_of_interest_by_name)

                of_interest_path = os.path.join(sub_dir_path, file_name.replace(SUFFIX, DOMAIN_OF_INTEREST + SUFFIX))
                with open(of_interest_path, 'w') as f:
                    # domain_of_interest_to_counts keys are exactly domains_of_interest_list but
                    # we need to keep the order in domains_of_interest_list
                    f.write('\n'.join(
                        domains_of_interest_to_names[domain] + '\t' + domain_of_interest_to_counts[domain] for domain in sorted_domain_of_interest_by_name))


                leftovers_path = os.path.join(sub_dir_path, file_name.replace(SUFFIX, LEFTOVERS + SUFFIX))
                sorted_leftovers = sorted(domain_leftovers_to_counts, key=domain_leftovers_to_counts.get, reverse=True)
                with open(leftovers_path, 'w') as f:
                    #print(domain_leftovers_to_counts)
                    '''for domain in domain_leftovers_to_counts:
                        print(domain+ '\t',)
                        print(domain_leftovers_to_counts[domain])
                        '''
                    print(leftovers_path)
                    f.write('\n'.join(
                        domain + '\t' +
                        domain_leftovers_to_counts[domain] for domain in sorted_leftovers))

                data_is_not_empty = True # to determine if the files should be written or not

                #copy_number file
                total_of_interest = sum(int(domain_of_interest_to_counts[domain]) for domain in domains_of_interest_list)
                copy_number_result += ','.join([sub_dir_name] + [domain_of_interest_to_counts[domain] for domain in domains_of_interest_list] + [str(total_of_interest)]) + '\n'

                #corrected_copy_number file
                total_of_interest_corrected = sum(int(domain_of_interest_to_corrected_counts[domain]) for domain in domains_of_interest_list)
                corrected_copy_number_result += ','.join([sub_dir_name] + [str(domain_of_interest_to_corrected_counts[domain]) for domain in domains_of_interest_list] + [str(total_of_interest_corrected)]) + '\n'

                #corrected_percent file
                corrected_percents = [float(domain_of_interest_to_corrected_counts[domain]) / total_of_interest_corrected for domain in domains_of_interest_list]
                corrected_percents_result += ','.join([sub_dir_name] + ['%.7f' % percent for percent in corrected_percents]) + '\n'

                #real_vs_leftovers file
                info_file_path = os.path.join(sub_dir_path, file_name.replace(COUNTS, INFO))
                with open(info_file_path) as f:
                    for line in f:
                        if TOTAL_SEQ_AFTER_FILTERS in line:
                            total_seq_after_filters = int(line.rstrip().split()[-1])
                            break
                total_leftovers = total_seq_after_filters - total_of_interest
                percent_of_interest = total_of_interest / float(total_seq_after_filters) * 100
                percent_leftovers = total_leftovers / float(total_seq_after_filters) * 100
                real_vs_leftovers_result += ','.join([of_interest_path, str(total_of_interest), '%.7f' % percent_of_interest, leftovers_path, str(total_leftovers), '%.7f' % percent_leftovers]) + '\n'

if data_is_not_empty:
    with open(copy_number_file, 'w') as f:
        f.write(copy_number_result)

    with open(corrected_copy_number_file, 'w') as f:
        f.write(corrected_copy_number_result)

    with open(corrected_percents_file, 'w') as f:
        f.write(corrected_percents_result)

    with open(real_vs_leftovers_file, 'w') as f:
        f.write(real_vs_leftovers_result)
else:
    logger.info('No relevant data for this run (no files will be written)')

logger.info('Done aggregating.')
end = time()
logger.info('Finished joining samples. Took {}'.format(measure_time(start, end)))