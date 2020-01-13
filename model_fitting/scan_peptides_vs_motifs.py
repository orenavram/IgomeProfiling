import datetime
import subprocess
import logging
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty


def calculate_pssm_thresholds(meme_path, cutoffs_path, faa_path, number_of_random_pssms, output_path, done_path,
                              argv='no_argv',
                              pssm_score_peptide='/groups/pupko/orenavr2/igomeProfilingPipeline/src/PSSM_score_Peptide/PSSM_score_Peptide'):

    if not os.path.exists(output_path):
        # TODO: any modules to load?
        cmd = f'{pssm_score_peptide} -pssm {meme_path} -pssm_cutoffs {cutoffs_path} -seq {faa_path} ' \
              f'-out {output_path} -NrandPSSM {number_of_random_pssms} -CalcPSSM_Pval'
        logger.info(f'{datetime.datetime.now()}: starting CalcPSSM_Pval. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)
    else:
        logger.info(f'{datetime.datetime.now()}: skipping scanning calculation as it is already exist at:\n{output_path}')

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(output_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_file_path', help='A path to a meme file with motifs against which a set of random peptides will be scanned')
    parser.add_argument('cutoffs_file_path', help='A path to a cutoffs file (peptide above cutoff? -> peptide is part of the motif')
    parser.add_argument('faa_file_path', help='A path to a faa file with peptides to scan against the pssms in the meme file')
    parser.add_argument('number_of_random_pssms', type=int, help='Number of pssm permutations')
    parser.add_argument('output_path', help='A path to which the Pvalues will be written to')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    calculate_pssm_thresholds(args.meme_file_path, args.cutoffs_file_path, args.faa_file_path,
                              args.number_of_random_pssms, args.output_path, args.done_file_path, argv=sys.argv)

