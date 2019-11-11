import gzip
import regex
import os
import datetime



def nnk_translate(dna_sequence, nnk_table):

    dna_sequence = dna_sequence
    seq_length = len(dna_sequence)
    aa_seq = ""
    for i in range(0, seq_length - 2, 3):
        codon = dna_sequence[i:i + 3]
        if codon in nnk_table:
            aa = nnk_table[codon]
        else:
            if codon == "TGA" or codon == "TAA":
                return 0 #"STOP_CODON"
            else:
                return 0 #"NOT_NNK"
        aa_seq += aa
    return aa_seq


def filter_reads(fastq_file, output_dir, barcode2samplename_path,
                 left_construct, right_construct, max_mismatches_allowed,
                 min_sequencing_quality, lib_types):

    nnk_table: {str: str} = {"CGT": "R", "CGG": "R", "AGG": "R",
                 "CTT": "L", "CTG": "L", "TTG": "L",
                 "TCT": "S", "TCG": "S", "AGT": "S",
                 "GCT": "A", "GCG": "A",
                 "GGT": "G", "GGG": "G",
                 "CCG": "P",
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

    logger.info(f'{datetime.datetime.now()}: Loading barcode2samplename file from:\n{barcode2samplename_path}')
    barcode_to_samplename = {}
    with open(barcode2samplename_path) as f:
        for line in f:
            if line.isspace():  # empty line
                continue
            barcode, sample_name = line.strip().split()
            if barcode in barcode_to_samplename:
                assert False, f'Barcode {barcode} belongs to more than one sample!!' # TODO: write to a global error log file
            barcode_to_samplename[barcode] = sample_name
    assert len(barcode_to_samplename) > 0, f'No barcodes were found in {barcode2samplename_path}'  # TODO: add informative error to log
    barcode_len = len(barcode)

    barcode2filenames = {}
    barcode2filehandlers = {}

    for barcode in barcode_to_samplename:
        sample_name = barcode_to_samplename[barcode] #f'{barcode}_{os.path.split(fastq_file)[-1]}'
        base_name = sample_name #f'{barcode}_{os.path.split(fastq_file)[-1]}'

        # create a sub dir for that barcode (under output_dir)
        subdir_path = f'{output_dir}/{sample_name}_{barcode}'
        os.makedirs(subdir_path, exist_ok=True)

        # add filenames
        barcode2filenames[barcode] = {}
        barcode2filenames[barcode]['fastq'] = f'{subdir_path}/{base_name}.fastq'  # translated seq # TODO: change to .gz
        barcode2filenames[barcode]['reads'] = f'{subdir_path}/{base_name}.fna'  # translated seq # TODO: change to .gz
        barcode2filenames[barcode]['peptides'] = f'{subdir_path}/{base_name}.faa'  # translated seq
        barcode2filenames[barcode]['filtered'] = f'{subdir_path}/{base_name}.filtration_log.txt'  # where in the filter seqs dropped

        # add filehandlers
        barcode2filehandlers[barcode] = {}
        barcode2filehandlers[barcode]['fastq'] = open(barcode2filenames[barcode]['fastq'], 'w')  # TODO: change to zip.open, 'wb'
        barcode2filehandlers[barcode]['reads'] = open(barcode2filenames[barcode]['reads'], 'w')  # TODO: change to zip.open, 'wb'
        barcode2filehandlers[barcode]['peptides'] = open(barcode2filenames[barcode]['peptides'], 'w')
        barcode2filehandlers[barcode]['filtered'] = open(barcode2filenames[barcode]['filtered'], 'w')

    log_f = open(f'{output_dir}/summary_log.txt', 'w')
    # counters:
    no_barcode = 0
    low_quality = 0
    too_many_mistakes = 0
    not_nnk = 0

    line_num = 0
    with gzip.open(fastq_file, 'rb') as fastq_f:
        for line in fastq_f:
            n_count = 0
            line_num += 1
            # ignore first line of each record
            seq = fastq_f.readline().decode("utf-8").rstrip()
            fastq_f.readline()  # ignore third line of each record
            quality = fastq_f.decode("utf-8").rstrip()

            barcode = seq[:barcode_len]
            if barcode not in barcode_to_samplename:
                no_barcode += 1
                continue

            barcode_quality = quality[:barcode_len]
            lowest_scored_character = min([ord(c) for c in barcode_quality])
            if lowest_scored_character < min_sequencing_quality:
                low_quality += 1
                barcode2filehandlers[barcode]["filtered"].write(f"{seq}\tlow quality barcode ({lowest_scored_character})\n")
                continue

            left_seq = seq[barcode_len: barcode_len + left_construct_length]
            if regex.match(left_seq + '{s<='f"{max_mismatches_allowed}"'}', left_construct) == "":
                too_many_mistakes += 1
                barcode2filehandlers[barcode]["filtered"].write(f"{seq}\tleft constract has more than {max_mismatches_allowed} mismatches\n")
                continue

            # TODO: make one file for short reads analysis and one for long ones
            right_seq = seq[barcode_len + left_construct_length:]
            match = regex.search(right_construct + '{s<='f"{max_mismatches_allowed}"'}', right_seq)
            if match != None:
                seq = right_seq[:match.span()[0]]
                n_count += seq.count('N')
                aa_seq = nnk_translate(seq, nnk_table)
                if aa_seq != 0 and n_count == 0:
                    if aa_seq.startswith('C') and aa_seq.endswith('C'):
                        libtype = f"C{len(aa_seq)-2}C"
                    else:
                        libtype = str(len(aa_seq))
                    barcode2filehandlers[barcode]["peptides"].write(f">LIB_{libtype}\n")
                else:
                    barcode2filehandlers[barcode]["filtered"].write(f"error aa_seq != 0 and n_count == 0\n")
            else:
                barcode2filehandlers[barcode]["filtered"].write(f"right_constract\n")



    for barcode in barcodes:
        barcode2filehandlers[barcode]["peptides"].close()
        barcode2filehandlers[barcode]["filtered"].close()

    log_f.close()


if __name__ == '__main__':
    from sys import argv
    print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_file', type=str, help='fastq file')
    parser.add_argument('output_dir', type=str, help='folder output')
    parser.add_argument('barcode2samplename', type=str, help='A path to the barcode to sample name file')
    # parser.add_argument('--gz', action='store_true', help='Turn on if the fastq is zipped (False by default)')
    parser.add_argument('--left_construct', type=str, default="AGGCGGCCAACGTGGC", help='SfilSite')
    parser.add_argument('--right_construct', type=str, default="GCCGCTGGGGCCGACC", help='RightConstract')
    parser.add_argument('--MistakeAllowed', type=int, default=1,
                        help='number of mistakes allowed in each construct')
    parser.add_argument('--min_sequencing_quality', type=int, default=38,
                        help='Minimum average sequencing threshold allowed after filtration'
                             'for more details, see: https://en.wikipedia.org/wiki/Phred_quality_score')
    parser.add_argument('--lib_types', type=str, default='C6C,6,C8C,8,C10C,10,C12C,12', help='CxC,x')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    filter_reads(args.fastq_file, args.output_dir, args.barcode2samplename,
                 args.left_construct, args.right_construct, args.MistakeAllowed,
                 args.min_sequencing_quality, args.lib_types)
