def peptides_intersection(filtered_reads_path, biological_condition, output_path, file_extension):
    peptides_files_counter = 0
    for path, dirs, files in os.walk(filtered_reads_path):
        for file in files:
            if biological_condition in path and file.endswith(file_extension):# and 'test' not in path:
                peptides_files_counter += 1
                file_path = os.path.join(path, file)
                logger.info(f'Loading {file_path}...')
                with open(file_path) as f:
                    motifs_frequency_counter: {int: int} = dict.fromkeys(range(1, 5), 0)
                    total_number_of_motifs = 0
                    for line in f:
                        samples_involved = set()
                        tokens = line.rstrip().split(',')
                        for token in tokens:
                            samples_involved.add(token.split('_')[1])
                        total_number_of_motifs += len(tokens)
                        motifs_frequency_counter[len(samples_involved)] += len(tokens)
                break
        if peptides_files_counter == 4:
            # relevant only for exp12. keep 5th sample away.
            break

    logger.info(f'peptides_files_counter={peptides_files_counter}')
    logger.info(f'total_number_of_motifs={total_number_of_motifs}')

    # save the results to a tab delimited file
    with open(output_path, 'a') as f:
        for frequency in sorted(motifs_frequency_counter):
            f.write(f'{frequency}\t')
            f.write(f'{motifs_frequency_counter[frequency]/total_number_of_motifs*100:.2f}\t')
            f.write(f'{biological_condition}\n')

import os
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('filtered_reads_path',
                        help='In this path each sample contains a folder in which there is a curated unique peptides file')
    parser.add_argument('biological_condition', help='the biological condition that its intersection should be quantified')
    parser.add_argument('output_path', help='a path to write the output to')
    parser.add_argument('--file_extension', default='cluster_to_combine.csv', help='the extension of the unique peptides file')
    args = parser.parse_args()


    os.makedirs(os.path.split(args.output_path)[0], exist_ok=True)
    with open(args.output_path, 'w') as f:
        f.write(f'intersection_size\tpercent\tmAb\n')

    for bc in args.biological_condition.split(','):
        logger.error(f'Parsing {bc}...')
        peptides_intersection(args.filtered_reads_path, bc, args.output_path, args.file_extension)


# python gershoni/src/IgOmeProfiling/motifs_intersection.py /Users/Oren/Dropbox/Projects/gershoni/IgOmeProfiling/exp12/pipeline_outputs/motif_inference_real b12,17b,21c,Herceptin /Users/Oren/Dropbox/Projects/gershoni/IgOmeProfiling/exp12/motifs_intersection.txt

