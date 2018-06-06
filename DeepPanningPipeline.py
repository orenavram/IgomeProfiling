from sys import argv


if len(argv) < 3:
    print('Usage: python '+argv[0]+' <fastq_path> <out_dir> <names_to_domains_of_interest_file> <lib_name>; <?mistakes_allowed> <?barcode_to_name (tab delimited)> <?code_to_barcode (tab delimited)>')
    # exit(-1)