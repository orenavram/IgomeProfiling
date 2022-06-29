import datetime
import subprocess
import logging
import os
import math
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, get_cluster_size_from_name, load_fasta_to_dict, \
    get_count_from, get_samples_from_name, get_unique_members_from_name


def get_clusters_sequences(motif_inference_output_path, biological_condition, sample_names,
                           cluster_names, cluster_rank, max_number_of_sequences_to_use):
    sample_paths = []
    for sample_name in sample_names:
        sample_paths.append(os.path.join(motif_inference_output_path, sample_name, 'unaligned_sequences'))

    unified_dict_sequence2header = {}  # unified_cluster
    for cluster_file_name in cluster_names:
        for sample_folder in sample_paths:
            if cluster_file_name in os.listdir(sample_folder):
                cluster_file_path = os.path.join(sample_folder, cluster_file_name)
                break
        else:
            raise ValueError(f'No cluster named {cluster_file_name} was found in the following dirs:\n' +
                             '\n'.join(sample_paths))
        sequence2header = load_fasta_to_dict(cluster_file_path, reverse=True)[0]  # don't need the other returned objects
        for sequence, header in sequence2header.items():
            if sequence in unified_dict_sequence2header:
                unified_header = unified_dict_sequence2header.pop(sequence)
                total_counts = get_count_from(header) + get_count_from(unified_header)
                header = unified_header[:unified_header.rindex('_')] + f'_{total_counts}'
            unified_dict_sequence2header[sequence] = header

    # TODO: should the counts be normalized by the number of samples involved? each sample contributes million reads
    # TODO: so 6 samples will contribute more than just 2 samples...
    # TODO: if so, we should divide HERE the total count (last token of the header) by the number of samples...
    unified_dict_header2sequence = {header: sequence for sequence, header in unified_dict_sequence2header.items()}

    result = ''
    for i, header in enumerate(sorted(unified_dict_header2sequence, key=get_count_from, reverse=True)):
        if i == max_number_of_sequences_to_use:
            break
        result += f'>{header}\n{unified_dict_header2sequence[header]}\n'

    number_of_unique_members = len(unified_dict_header2sequence)
    cluster_size = sum(get_count_from(header) for header in unified_dict_header2sequence)
    result_file_name = f'{biological_condition}_clusterRank_' \
                       f'{str(cluster_rank).zfill(4)}_uniqueMembers_' \
                       f'{"top" if number_of_unique_members >= max_number_of_sequences_to_use else ""}' \
                       f'{min(number_of_unique_members, max_number_of_sequences_to_use)}_' \
                       f'clusterSize_{cluster_size:.2f}.faa'
    return result, result_file_name


def remove_consensus_name(clusters):
    #remove consensus sequence so we have the exact cluster (file) name
    return [cluster[cluster.index('_')+1:] for cluster in clusters]


def get_samples_number_build_cluster(clusters):
    return len(set([get_samples_from_name(cluster) for cluster in clusters]))


def sort_clusters(clusters):
    sort_by_num_samples = get_samples_number_build_cluster(clusters)
    sort_by_unique_memebers = sum(get_unique_members_from_name(cluster) for cluster in clusters)*(-1)
    sort_by_cluster_size = sum(get_cluster_size_from_name(cluster) for cluster in clusters)
    return (sort_by_num_samples, sort_by_cluster_size , sort_by_unique_memebers)

def unite_clusters(motif_inference_output_path, meme_file, biological_condition, sample_names,
                   max_number_of_members_per_cluster, output_path, done_path, aln_cutoff, pcc_cutoff, 
                   sort_cluster_to_combine_only_by_cluster_size, min_number_samples_build_cluster,
                   unite_pssm_script_path='./UnitePSSMs/UnitePSSMs', argv='no_argv'):

    clusters_to_combine_path = os.path.join(output_path, 'cluster_to_combine.csv')
    if not os.path.exists(clusters_to_combine_path):
        # TODO: any modules to load?
        cmd = f'{unite_pssm_script_path} -pssm {meme_file} -out {clusters_to_combine_path} ' \
              f'-aln_cutoff {aln_cutoff} -pcc_cutoff {pcc_cutoff}'
        logger.info(f'{datetime.datetime.now()}: starting UnitePSSMs. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(clusters_to_combine_path)

    logger.info(f'Result file is at {clusters_to_combine_path}')
    clusters_to_combine = []
    with open(clusters_to_combine_path) as f:
        for line in f:
            cluster_names = line.rstrip().split(',')
            clusters_to_combine.append(cluster_names)

    logger.info(f'Sorting clusters by rank...')
    # sort the sublist such that the first one will contain the highest copy number, etc...
    if sort_cluster_to_combine_only_by_cluster_size:
            clusters_to_combine.sort(key=lambda clusters: sum(get_cluster_size_from_name(cluster) for cluster in clusters), reverse=True)
    else:
        clusters_to_combine.sort(key=lambda clusters: sort_clusters(clusters), reverse=True)
    sorted_clusters_to_combine_path = clusters_to_combine_path.replace('cluster_to_combine', 'sorted_cluster_to_combine')
    with open(sorted_clusters_to_combine_path, 'w') as f:
        for cluster_names in clusters_to_combine:
            num_of_samples_per_cluster = get_samples_number_build_cluster(cluster_names) 
            if num_of_samples_per_cluster >= min_number_samples_build_cluster:
                f.write(','.join(cluster_names)+'\n')

    unaligned_sequences_path = os.path.join(output_path, 'unaligned_sequences')
    os.makedirs(unaligned_sequences_path, exist_ok=True)

    for cluster_rank in range(len(clusters_to_combine)):
        if cluster_rank % 25 == 0:
            logger.info(f'Merging sequences of the cluster ranked {cluster_rank}')

        clusters_sequences, cluster_file_name = get_clusters_sequences(motif_inference_output_path, biological_condition,
                                                                       sample_names, remove_consensus_name(clusters_to_combine[cluster_rank]),
                                                                       cluster_rank, max_number_of_members_per_cluster)
        with open(os.path.join(unaligned_sequences_path, cluster_file_name), 'w') as f:
            f.write(clusters_sequences)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('motif_inference_output_path',
                        help='A path to a folder in which each subfolder contains a different sample analysis')
    parser.add_argument('meme_file_path', help='A path to a meme file')
    parser.add_argument('biological_condition', help='A biological condition that its motifs will be unified')
    parser.add_argument('sample_names', help='Sample names to apply over the "motif unification". More than one '
                                                 'sample name should be separated by commas but no spaces. '
                                                 'For example: 17b_01,17b_03,17b_05')
    parser.add_argument('max_number_of_members_per_cluster', type=int,
                        help='How many members (at most) should be taken to each cluster')
    parser.add_argument('output_path', help='A path in which a new subfolder with the united motifs will be written to')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--aln_cutoff', default='24', help='The cutoff for pairwise alignment score to unite motifs of BC')
    parser.add_argument('--pcc_cutoff', default='0.7', help='Minimal PCC R to unite motifs of BC')
    parser.add_argument('--sort_cluster_to_combine_only_by_cluster_size', action='store_true', help='Sort the clusters only by the cluster size')
    parser.add_argument('--min_number_samples_build_cluster_per_BC', type=str, default=1, help='Keep only clusters that build from X minimum number of samples')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    unite_clusters(args.motif_inference_output_path, args.meme_file_path, args.biological_condition,
                   args.sample_names.split(','), args.max_number_of_members_per_cluster, args.output_path, args.done_file_path,
                   args.aln_cutoff, args.pcc_cutoff, args.sort_cluster_to_combine_only_by_cluster_size,
                   args.min_number_samples_build_cluster, argv=sys.argv)
