def load_sequences(file_path):
    seqs = []
    with open(file_path) as f:
        for line in f:
            if not line.startswith('>'):
                seqs.append(line.rstrip())
    logger.info(f'Number of peptides in {os.path.split(file_path)[-1]} is {len(seqs)}')
    return seqs

def peptides_intersection(filtered_reads_path, biological_condition, output_path, file_extension):
    seqs_frequency_counter: {str: int} = {}
    peptides_files_counter = 0
    for path, dirs, files in os.walk(filtered_reads_path):
        for file in sorted(files):
            if biological_condition in path and file.endswith(file_extension):
                peptides_files_counter += 1
                file_path = os.path.join(path,file)
                logger.info(f'Loading {file_path}...')
                for seq in load_sequences(file_path):
                    seqs_frequency_counter[seq] = seqs_frequency_counter.get(seq, 0) + 1
                break
        if peptides_files_counter == 4:
            # relevant only for exp12. keep 5th sample away.
            break

    logger.info(f'peptides_files_counter={peptides_files_counter}')
    # frequency counter of the values
    intersection_size_to_count: {int: int} = dict.fromkeys(range(1, peptides_files_counter+1), 0)
    # foreach 1<=k<=4 a mapping from k to the set of sequences that were shared between k samples
    shared_in_k_samples: {int: str} = {i: set() for i in range(1, 5)}
    for seq in seqs_frequency_counter:
        intersection_size = seqs_frequency_counter[seq]
        intersection_size_to_count[intersection_size] += 1
        shared_in_k_samples[intersection_size].add(seq)  # TODO: what to do with this?

    with open(f'{os.path.split(output_path)[0]}/unshared_peptides_{biological_condition}.txt', 'w') as f:
        f.write('\n'.join(shared_in_k_samples[1]))

    total_number_of_peptides = sum(intersection_size_to_count.values())
    assert sum([key*intersection_size_to_count[key] for key in intersection_size_to_count]) == sum(seqs_frequency_counter.values())
    assert peptides_files_counter > 0

    # save the results to a tab delimited file
    with open(output_path, 'a') as f:
        for intersection_size in sorted(intersection_size_to_count):
            f.write(f'{intersection_size}\t')
            f.write(f'{intersection_size_to_count[intersection_size]/total_number_of_peptides*100}\t')
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
    parser.add_argument('--file_extension', default='rpm.faa', help='the extension of the unique peptides file')
    args = parser.parse_args()

    os.makedirs(os.path.split(args.output_path)[0], exist_ok=True)
    with open(args.output_path, 'w') as f:
        f.write(f'intersection_size\tpercent\tmAb\n')

    for bc in args.biological_condition.split(','):
        logger.error(f'Parsing {bc}...')
        peptides_intersection(args.filtered_reads_path, bc, args.output_path, args.file_extension)


# python gershoni/src/IgOmeProfiling/peptides_intersection.py /Users/Oren/Dropbox/Projects/gershoni/IgOmeProfiling/exp12/pipeline_outputs/first_phase_output_real b12,17b,21c,Herceptin /Users/Oren/Dropbox/Projects/gershoni/IgOmeProfiling/exp12/peptides_intersection.txt

