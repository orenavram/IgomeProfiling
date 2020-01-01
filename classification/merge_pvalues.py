import logging
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty


def get_consensus_sequences_from_meme(meme_path):
    result = []
    with open(meme_path) as f:
        for line in f:
            if not line.startswith('MOTIF'):
                continue
            # "MOTIF CLKGASFLAC_17b_clusterRank_0_uniqueMembers_339_clusterSize_2659173.71.faa"
            result.append(line.split()[1].split('_')[0])

    return result


def get_results(path):
    pvalues = []
    hits = []
    with open(path) as f:
        for line in f:
            if line.startswith('##'):
                # "## PSSM_name	p_Value	True_Hits: num_of_hits"
                continue
            # "CLKGASFLAC_17b_clusterRank_0_uniqueMembers_339_clusterSize_2659173.71.faa\t0.01\tTrue_Hits: 118"
            line_tokens = line.split('\t')
            # consensus.append(line_tokens[0].split('_')[0])
            pvalues.append(line_tokens[1])
            hits.append(line_tokens[2].split()[-1])

    return pvalues, hits


def aggregate_pvalues_results(meme_path, scanning_results_dir_path, aggregated_pvalues_path,
                              aggregated_hits_path, done_path, argv='no argv'):

    all_consensuses = get_consensus_sequences_from_meme(meme_path)
    pvalues_f = open(aggregated_pvalues_path, 'w')
    hits_f = open(aggregated_hits_path, 'w')

    #header
    pvalues_result = hits_result = f'sample_name,{",".join(all_consensuses)}\n'
    # consensus2pvalues_column = {consensus:[] for consensus in all_consensuses}
    # consensus2hits_column = {consensus:[] for consensus in all_consensuses}
    for file_name in sorted(os.listdir(scanning_results_dir_path)):
        if file_name.endswith('100.txt'):
            raise TypeError

        if file_name.endswith('00.txt'):
            # next sample is starting
            pvalues_f.write(f'{pvalues_result.rstrip().rstrip(",")}\n')
            hits_f.write(f'{hits_result.rstrip().rstrip(",")}\n')
            sample_name = file_name.split('_peptides')[0]
            pvalues_result = hits_result = f'{sample_name},'

        pvalues, hits = get_results(os.path.join(scanning_results_dir_path, file_name))
        print(len(pvalues))
        pvalues_result += ','.join(pvalues) + ','
        hits_result += ','.join(hits) + ','

    pvalues_f.write(f'{pvalues_result.rstrip(",")}\n')
    hits_f.write(f'{hits_result.rstrip(",")}\n')

    pvalues_f.close()
    hits_f.close()

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(aggregated_pvalues_path)
    verify_file_is_not_empty(aggregated_hits_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_path', help='A path to a (MEME) file with the motifs against which the random peptides were scanned (in silico)')
    parser.add_argument('scanning_results_dir_path', help='A path in which each file contains a hits/pvals computation. '
                                                         'These files will be aggregated into one table.')
    parser.add_argument('aggregated_pvalues_path', help='A path to which the Pvalues table will be written to')
    parser.add_argument('aggregated_hits_path', help='A path to which the hits table will be written to')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script was finished running successfully.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    aggregate_pvalues_results(args.meme_path, args.scanning_results_dir_path, args.aggregated_pvalues_path,
                              args.aggregated_hits_path, args.done_file_path, argv=sys.argv)






'''
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('inputfolder', type=str, help='input folder')
parser.add_argument('outputfolder', type=str, help='folder output')
parser.add_argument('sample_to_biological', type=str, help='sample_to_biological file')
args = parser.parse_args()

inputfolder = args.inputfolder.rstrip('/') + '/'
outputfolder = args.outputfolder.rstrip('/') + '/'
sample_to_biological = args.sample_to_biological

biological_condition_of_folder = inputfolder.rstrip('/')
biological_condition_of_folder = biological_condition_of_folder[biological_condition_of_folder.rfind('/') + 1:].lower()
sample2biologicalcondition = {}

with open(sample_to_biological, 'r') as f:
    for line in f:
        sample, bio = line.rstrip().split("\t")
        sample2biologicalcondition[sample.lower().replace('_', '-')] = bio.lower()


def get_label(sample):
    if biological_condition_of_folder in sample2biologicalcondition[sample.lower()]:  # and "05" not in samples:
        return biological_condition_of_folder
    else:
        return "other"


motif_to_pvalues = {}
motif_to_hits = {}
samples = []
sample = ""
for file in sorted(os.listdir(inputfolder)):
    if not file.endswith(
            'S100.txt'):  # or "-05" in file: # need to change the -05 its for exp12 when 05 is the test sample
        continue
    current_sample = file[file.find("sample_") + len("sample_"):file.rfind('_')]
    if sample != current_sample:
        sample = current_sample
        samples.append(sample)
    with open(inputfolder + file, 'r') as f:
        for line in f:
            if line.startswith("##") or line.strip() == "":
                # skip header
                continue
            line = line.rstrip().split()
            motif = line[0][:line[0].find('_')]
            motif_to_pvalues[motif] = motif_to_pvalues.get(motif, [])
            motif_to_hits[motif] = motif_to_hits.get(motif, [])
            motif_to_pvalues[motif].append(float(line[1]))
            motif_to_hits[motif].append(float(line[3]) + 1)

check_for_keys = []
for i in range(len(samples)):
    if biological_condition_of_folder in sample2biologicalcondition[samples[i].lower()]:  # and "05" not in samples[i]:
        check_for_keys.append(i)

key_to_del = []
for motif in motif_to_pvalues:
    j = 0
    for i in check_for_keys:
        pvalues = motif_to_pvalues[motif][i]
        if pvalues > 0.05:
            j += 1
    if j == len(check_for_keys):
        key_to_del.append(motif)

for key in key_to_del:
    del motif_to_pvalues[key]

with open(outputfolder + f"hits/{biological_condition_of_folder}_hits.csv", 'w') as f:
    f.write(f"{','.join(motif_to_hits.keys())},label,sample_name\n")
    for i in range(len(samples)):
        for key in motif_to_hits:
            f.write(f"{motif_to_hits[key][i]},")
        f.write(f"{get_label(samples[i])},{samples[i]}\n")

with open(outputfolder + f"pvalues/{biological_condition_of_folder}_pvalues.csv", 'w') as f:
    f.write(f"{','.join(motif_to_pvalues.keys())},label,sample_name\n")
    for i in range(len(samples)):
        for key in motif_to_pvalues:
            f.write(f"{motif_to_pvalues[key][i]},")
        f.write(f"{get_label(samples[i])},{samples[i]}\n")
'''
