import sys
import csv  
import os
import math
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import load_table_to_dict


def get_rpm_factor(parsed_fastq_results, sample_name_of_barcode, delimiter='\t'):
    path_rpm_factor = f'{parsed_fastq_results}/rpm_factors.txt'
    file_rpm_factor = open(path_rpm_factor, 'r')
    for line in file_rpm_factor.readlines():
        if line.isspace():  # empty line
            continue
        sample_name, count_peptides, rpm = line.strip().split(delimiter)
        if sample_name == sample_name_of_barcode:
            return str(rpm)


def get_number_unique_peptides(parsed_fastq_results, sample_name):
    path_unique_rpm_file = f'{parsed_fastq_results}/{sample_name}/{sample_name}_unique_rpm.faa'
    num_unique_peptides = str(math.ceil(len(open(path_unique_rpm_file, 'r').readlines()) / 2))
    return num_unique_peptides


def get_num_peptides_to_barcode(lines, barcode):
    for line in lines:
        if line.startswith(barcode):
            return line.split()[-1]


def get_count_peptides_before_and_after_filtration(parsed_fastq_results, barcode, num_barcodes):
    path_summary_log_txt = f'{parsed_fastq_results}/summary_log.txt'
    file_summary_log_txt = open(path_summary_log_txt, 'r')
    lines = file_summary_log_txt.readlines()
    before_filtration = ''
    after_filtration = ''
    for ind,line in enumerate(lines):
        if line.startswith('Total number of dna sequences with a LEGAL barcode'):
           before_filtration = get_num_peptides_to_barcode(lines[ind + 1 : ind + 1 + num_barcodes], barcode)
        elif line.startswith('Translated peptides (per barcode)'):
            after_filtration = get_num_peptides_to_barcode(lines[ind + 1 : ind + 1 + num_barcodes], barcode)
            break
        else:
            continue
    return [before_filtration, after_filtration]
    

def summary_reads(barcode2samplename, parsed_fastq_results, done_file_path, name_summary_file_reads, no_calculate_rpm, argv):
    barcode2samplename_dict = load_table_to_dict(barcode2samplename, 'Barcode {} belongs to more than one sample_name!!')
    
    header = ["barcode", "sample name", "Before filtration per barcode", "After filtration per barcode", "Unique peptides per sample", "1 RPM per sample"]
    with open(f'{parsed_fastq_results}/{name_summary_file_reads}', 'w') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)
        # write the data
        for barcode in barcode2samplename_dict:
            data = []
            sample_name = barcode2samplename_dict[barcode]
            # Bracode
            data.append(barcode)

            # Sample name
            data.append(sample_name)

            # Before filtration and After filtration
            data = data + get_count_peptides_before_and_after_filtration(parsed_fastq_results, barcode, len(barcode2samplename_dict))
            
            # Unique peptides
            data.append(get_number_unique_peptides(parsed_fastq_results, sample_name))
            
            # 1 RPM
            rpm_factor = '' if no_calculate_rpm else get_rpm_factor(parsed_fastq_results, sample_name)
            data.append(rpm_factor)

            writer.writerow(data)
    
    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}', flush=True)
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('barcode2samplename', type=str, help='A path to the barcode to sample name file')
    parser.add_argument('parsed_fastq_results', type=str, help='output folder')
    parser.add_argument('done_file_path', type=str, help='A path to a file that signals that the summary reads finished running successfully')
    parser.add_argument('--name_summary_file_reads', default='summary_log_reads.csv', type=str, help='A name for summary file summary all reads.')
    parser.add_argument('--no_calculate_rpm', action='store_true', help='Disable normalize counts to "reads per million" (sequence proportion x 1,000,000)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    summary_reads(args.barcode2samplename, args.parsed_fastq_results, args.done_file_path, args.name_summary_file_reads, args.no_calculate_rpm, sys.argv)
