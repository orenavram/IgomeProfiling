from sys import argv


if len(argv) < 3:
    print('Usage: python '+argv[0]+' <fastq_path> <out_dir> <domain_name_to_domain_sequence_file> <lib_name>; <?mistakes_allowed> <?sample_barcode_to_sample_name (tab delimited)> <?code_to_barcode (tab delimited)>')
    # exit(-1)