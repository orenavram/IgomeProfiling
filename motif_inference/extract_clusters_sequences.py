import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, get_count_from


def load_member_prefix_to_record_dict(fasta_file, prefix_length=10):
    # Why 10 by default? see load_clusters_to_members_dict documentation
    member_prefix_to_record = {}  # short header to seq
    with open(fasta_file) as f:
        for header in f:
            sequence = f.readline()
            member_prefix_to_record[header[:prefix_length]] = f'{header}{sequence}'

    return member_prefix_to_record


def load_clusters_to_members_dict(clstr_file, member_prefix_to_record, prefix_length=10):
    cluster_to_members_records = {}  # short header to header
    with open(clstr_file) as f:
        for line in f:
            if line.startswith('>'):
                # clstr file cluster headers looks like this:
                # ">Cluster $i"
                cluster = line.split()[1]
                cluster_to_members_records[cluster] = []
            else:
                # clstr file member looks like this:
                # "0	12aa, >seq_7_lib_12_len_12... *"
                member = line.split()[2]
                # CD-hit uses only 20 chars for sequnce header in the clstr file so the header is usually trimmed (and
                # there is a "..." afterward). Thus, it's worthless to look for for the whole header in the fasta file
                # when the headers are longer than 20 chars. Just look for its prefix!
                # 10 is (hopefully) long enough to be unique and short enough to stop before the "..."
                member_prefix = member[:prefix_length]  # min(prefix_length, line.index('...'))
                cluster_to_members_records[cluster].append(member_prefix_to_record[member_prefix])

    return cluster_to_members_records


def extract_sequence_counts_from_record(rec):
    # Example record:
    # >someHeaderInWhichTheCountsReside\nHGKTGASFLQ
    header = rec.split('\n')[0]
    return get_count_from(header)


def extract_cluster_size_from_records(records_list):
    return sum(extract_sequence_counts_from_record(rec) for rec in records_list)


def extract_clusters_sequences(fasta_file, clstr_file, output_dir, done_path,
                               max_num_of_sequences_to_keep, file_prefix, argv='no argv'):

    verify_file_is_not_empty(fasta_file)
    verify_file_is_not_empty(clstr_file)

    os.makedirs(output_dir, exist_ok=True)

    member_prefix_to_record = load_member_prefix_to_record_dict(fasta_file)

    cluster_to_members_records = load_clusters_to_members_dict(clstr_file, member_prefix_to_record)

    logger.info(f'{datetime.datetime.now()}: Writing clusters sequences...')

    trimmed_clusters = set()
    # sort records of each cluster by their size and keep only first $max_num_of_sequences_to_keep records
    for cluster in cluster_to_members_records:
        # sort cluster members by their "strength", i.e., counts
        cluster_to_members_records[cluster].sort(key=extract_sequence_counts_from_record, reverse=True)
        if len(cluster_to_members_records[cluster])>max_num_of_sequences_to_keep:
            # discard (in-place) all sequences above the maximum required number
            cluster_to_members_records[cluster][max_num_of_sequences_to_keep:] = []
            trimmed_clusters.add(cluster)

    max_number_of_leading_zeros = len(str(len(cluster_to_members_records)))

    sorted_clusters_by_size = sorted(cluster_to_members_records, reverse=True,
                                     key=lambda cluster: extract_cluster_size_from_records(cluster_to_members_records[cluster]))

    if file_prefix != '':
        file_prefix += '_'

    for i, cluster in enumerate(sorted_clusters_by_size):

        cluster_rank = str(i).zfill(max_number_of_leading_zeros)
        number_of_unique_members = min(len(cluster_to_members_records[cluster]), max_num_of_sequences_to_keep)
        cluster_counts = extract_cluster_size_from_records(cluster_to_members_records[cluster])

        filename = f'{file_prefix}clusterRank_' \
                   f'{cluster_rank}_uniqueMembers_' \
                   f'{"top" if cluster in trimmed_clusters else ""}' \
                   f'{number_of_unique_members}_' \
                   f'clusterSize_{cluster_counts:.2f}.faa'  # take only 2 digits after the floating point

        with open(os.path.join(output_dir, filename), 'w') as f:
            f.write(''.join(record for record in cluster_to_members_records[cluster]))

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help='A fasta file to collapse for unique sequences and their counts')
    parser.add_argument('clstr_file', help='A .clstr file that details the cluster in the fasta_file')
    parser.add_argument('output_dir', help='A folder in which each cluster will be written as a separate file.')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script was finished running successfully.')
    parser.add_argument('--max_num_of_sequences_to_keep', type=int, default=1000,
                        help='How many sequences (at most) to include in each file? (extra sequences will not be '
                             'extracted. A very high number might not extract anything)')
    parser.add_argument('--file_prefix', default='',
                        help='An optional prefix str to add to each of the extracted clusters. For example, sample name')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging

    if args.verbose:
      logging.basicConfig(level=logging.DEBUG)
    else:
      logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    extract_clusters_sequences(args.fasta_file, args.clstr_file, args.output_dir, args.done_file_path,
                               args.max_num_of_sequences_to_keep, args.file_prefix, sys.argv)
