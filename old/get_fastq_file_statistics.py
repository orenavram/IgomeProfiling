def get_stats(fastq_path, output_path_prefix, upstream_construct,
              downstream_construct, mistake_allowed, gz):
    import regex
    import gzip

    nnk_translation = {'CGT': 'R', 'CGG': 'R', 'AGG': 'R',
                       'CTT': 'L', 'CTG': 'L', 'TTG': 'L',
                       'TCT': 'S', 'TCG': 'S', 'AGT': 'S',
                        'GCT': 'A', 'GCG': 'A',
                        'GGT': 'G', 'GGG': 'G',
                        'CCT': 'P', 'CCG': 'P',
                        'ACT': 'T', 'ACG': 'T',
                        'CAG': 'Q',
                        'TAG': 'q',
                        'GTT': 'V', 'GTG': 'V',
                        'AAT': 'N',
                        'GAT': 'D',
                        'TGT': 'C',
                        'GAG': 'E',
                        'CAT': 'H',
                        'ATT': 'I',
                        'AAG': 'K',
                        'ATG': 'M',
                        'TTT': 'F',
                        'TGG': 'W',
                        'TAT': 'Y'}
    libraries = ['6', '8', '10', '12', 'C6C', 'C8C', 'C10C', 'C12C', '14']

    # count for each random dna fragment its occurrences (for each library)
    library_to_dna_to_count = {library: {} for library in libraries}

    # how many random dna sequences appear k times for each k (for each library)
    library_to_repetitions_frequency_counter = {library: {} for library in libraries}

    open_function = gzip.open if gz else open
    with open_function(fastq_path, 'r' + 'b'*gz) as f:
        f.readline()
        sample_barcode = f.readline()[:8]

    f_filtered_out = open(os.path.join(output_path_prefix, f'filtered_out.txt'), 'w')
    with open_function(fastq_path, 'r' + 'b'*gz) as f:
        for line in f:
            # skip first line, read the second one, in which the read sequence appears
            line = f.readline()
            if gz:
                line = line.decode("utf-8")
            line = line[len(sample_barcode):].rstrip() # remove barcode and trailing \n

            f.readline()  # skip third line
            f.readline()  # skip fourth line
            match = regex.match(f'({upstream_construct})' + '{s<=' + f'{mistake_allowed}' + '}' +
                                '(([ACGT][ACGT][GT]){6,})' +  # this is an NNK format
                                f'({downstream_construct})' + '{s<=' +f'{mistake_allowed}' +'}', line)
            # An example regex is '(AGGCGGCCAACGTGGC){s<=1}(.+)(GCCGCTGGGGCCGACC){s<=1}'

            if not match:
                if not regex.match(f'({upstream_construct})' + '{s<=' + f'{mistake_allowed}' + '}', line):
                    f_filtered_out.write(f'Too many mismatches in UP stream construct\t{upstream_construct}\n')
                elif not regex.search(f'({downstream_construct})' + '{s<=' +f'{mistake_allowed}' +'}', line):
                    f_filtered_out.write(f'Too many mismatches in DOWN stream construct\t{downstream_construct}\n')
                else:
                    random_dna = regex.match(f'({upstream_construct})' + '{s<=' + f'{mistake_allowed}' + '}' +
                                '(.+)' + f'({downstream_construct})' + '{s<=' +f'{mistake_allowed}' +'}', line).group(2)
                    f_filtered_out.write(f'Random dna is not in legal NNK library format (length={len(random_dna)})\t{random_dna}\n')
                continue

            random_dna = match.group(2) # extract the random sequence between the upstream and downstream constructs

            # if len(random_dna)%3 != 0:
            #     logger.error(f'SKIPPING... Random dna is not in NNK format (cannot be divided by 3):\n{random_dna}')
            #     continue

            # try:
            random_aa = ''.join(nnk_translation[random_dna[i:i+3]] for i in range(0, len(random_dna), 3))
            # except:
            #     logger.info(f'SKIPPING... Random dna is not in NNK format (K is A or C). This is its 3, 6, 9 etc dna base pairs:\n{random_dna[2::3]}')
            #     continue
            aa_len = len(random_dna)//3

            try:
                if random_aa.startswith('C') and random_aa.endswith('C') and aa_len > 6:  # "C4C" should be saved as "6"
                    # Cys loop construct
                    library_to_dna_to_count[f'C{aa_len-2}C'][random_dna] = library_to_dna_to_count[f'C{aa_len-2}C'].get(random_dna, 0) + 1
                else:
                    library_to_dna_to_count[f'{aa_len}'][random_dna] = library_to_dna_to_count[f'{aa_len}'].get(random_dna, 0) + 1
            except:
                f_filtered_out.write(f'Illegal random library of length {len(random_dna)}\t{random_aa}\t{random_dna}\n')
                continue

    f_filtered_out.close()

    for lib in library_to_dna_to_count:
        for dna in library_to_dna_to_count[lib]:
            k = library_to_dna_to_count[lib][dna]
            library_to_repetitions_frequency_counter[lib][k] = library_to_repetitions_frequency_counter[lib].get(k, 0) + 1

    for lib in library_to_dna_to_count:
        with open(os.path.join(output_path_prefix, f'{lib}_counts.csv'), 'w') as f:
            total_reads = sum(library_to_dna_to_count[lib].values())
            f.write('random_dna\tcounts\tpercents\trandom_peptide\n')
            for dna in sorted(library_to_dna_to_count[lib], key=library_to_dna_to_count[lib].get, reverse=True):
                f.write(f'{dna}\t'
                        f'{library_to_dna_to_count[lib][dna]}\t'
                        f'{library_to_dna_to_count[lib][dna]/total_reads*100:.6f}%\t'
                        f'{"".join(nnk_translation[dna[i:i+3]] for i in range(0, len(dna), 3))}\n')

        with open(os.path.join(output_path_prefix, f'{lib}_lengths_distribution.csv'), 'w') as f:
            f.write('k\tnumber of reads that appear k times\tpercents\n')
            total_reads = sum(library_to_repetitions_frequency_counter[lib].values())
            for k in sorted(library_to_repetitions_frequency_counter[lib]):
                f.write(f'{k}\t'
                        f'{library_to_repetitions_frequency_counter[lib][k]}\t'
                        f'{library_to_repetitions_frequency_counter[lib][k]/total_reads*100:.6f}%\n')








if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import os

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('first_phase_output_dir', help='path to a first_phase_output_dir')
        # parser.add_argument('output_path_prefix', help='a prefix to a path in which the outputs will be written')
        # parser.add_argument('sample_barcode', help='sample barcode')
        parser.add_argument('--upstream_construct', help='upstream sequence construct', default='AGGCGGCCAACGTGGC')
        parser.add_argument('--downstream_construct', help='downstream sequence construct', default='GCCGCTGGGGCCGACC')
        parser.add_argument('--mistake_allowed', help='max num of mismatches tolerance', default=1, type=int)
        parser.add_argument('--output_path', help='a folder to which the output should be written', default='', type=int)
        parser.add_argument('-gz', '--gz', help='Increase output verbosity', action='store_true')
        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        import logging
        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger('main')

        nuc_fastq_suffix = '..NUC.fastq'
        if args.gz:
            nuc_fastq_suffix += '.gz'
        len_suffix = len(nuc_fastq_suffix)

        for path, dirs, files in os.walk(args.first_phase_output_dir):
            for file in files:
                if file.endswith(nuc_fastq_suffix):
                    fastq_path = os.path.join(path, file)
                    logger.info(f'Handling fastq at:\n{fastq_path}')
                    sample_name = os.path.split(path.rstrip('/'))[-1]
                    if args.output_path:
                        output_path = os.path.join(args.output_path, f'{sample_name}_stats')
                    else:
                        output_path = os.path.join(path, f'{sample_name}_stats')
                    os.makedirs(output_path, exist_ok=True)
                    get_stats(fastq_path, output_path, args.upstream_construct,
                              args.downstream_construct, args.mistake_allowed, args.gz)
