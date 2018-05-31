import gzip
from sys import argv
import os
import regex

# CONSTANTS:
fth1_anealed_site = "AAGTAGGGGATCCAGG"  # mistakes allowed
fth1_anealed_antisense = "TCTAGAGCCGACCGCGA"  # no mistakes allowed

sample_barcode_length = 8
fth1_anealed_site_length = len(fth1_anealed_site)
VarLength = 12;
fth1_anealed_antisense_length = len(fth1_anealed_antisense)
var_starts = sample_barcode_length + fth1_anealed_site_length
antisense_starts = var_starts + VarLength

FILTERED = 'filtered'
INFO = 'info'
FASTQ = 'fastq'
FS = 'fs'
COUNTS = 'counts'

TOTAL_SEQ_AFTER_FILTERS = 'TOTAL SEQ AFTER FILTERS'
DOMAIN_OF_INTEREST = '.domain_of_interest'
LEFTOVERS = '.domain_leftovers'
suffix = '.txt'

TooManyMistakes = 'TOO_MANY_MISTAKES'#'TooManyMistakes'
Seq_with_errors_in_fth1_annealed_anti_sense = 'Seq_with_errors_in_fth1_annealed_anti_sense'
Wrong_Seq_With_N = 'Wrong_Seq_With_N'
TotalSeq = 'TotalSeq'
TotalSeqsOK = 'TotalSeqsOK'

def load_table_to_dict(file_to_read, d = None, delimiter ='\t', backwards = 0, open_operator=open):
    '''
    parse a delimited file with two columns to a dictionary
    if backwards == 0 keys are taken from the first column o.w., from the second
    open_operator can also be set to gzip.open if the file is zipped
    '''
    if d == None:
        d = {}
    with open_operator(file_to_read) as f:
        for line in f:
            if line != '' and not line.isspace():
                item1, item2 = line.rstrip().split(delimiter)
                d[item1] = item2
    if backwards:
        d = dict((d[key],key) for key in d)
    return d #, dict([(barcode_to_name[barcode], barcode) for barcode in barcode_to_name])

def diff(s1, s2):
    assert len(s1) == len(s2)
    return sum([s1[i]!=s2[i] for i in range(len(s1))])

def load_column_to_list(file_path):
    with open(file_path) as f:
        return [line.rstrip() for line in f]

def measure_time(start, end):
    hours = int(end - start) / 3600
    minutes = (int(end - start) % 3600) / 60
    seconds = int(end - start) % 60
    return '%d:%d:%d hours.' % (hours, minutes, seconds)

