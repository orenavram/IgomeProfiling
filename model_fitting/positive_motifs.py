import logging
import os
import sys
from matplotlib.pyplot import axis
import pandas as pd
import numpy as np
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)


def normalize_max(df, max_value):
    return (df / max_value).replace(np.nan, 0)


def normalize_max_min(df, max_value, min_value):
    return ((df - min_value) / (max_value - min_value)).replace(np.nan, 0)


def normalize_log(df, rank_method):
    if rank_method == 'hits':
        return np.log2(df + 1)
    if rank_method == 'shuffles':
        return -np.log2(1- df + 0.01)
    return df


def min_max_fixed(df, max_value, min_value):
    df[df > max_value] = max_value
    df[df < min_value] = min_value
    return df


def is_artifact(df, bio_cond, invalid_mix):
    artifact_motifs = []
    motifs = list(df.columns)
    motifs.remove('sample_name')
    motifs.remove('label')
    for motif in motifs:
        bc_min = df.loc[df['label'] == bio_cond, motif].min()
        mixed_samples = list(df.loc[(df['label'] == 'other') & (df[motif] >= bc_min), 'sample_name'])
        is_artifact = any(invalid_mix in s for s in mixed_samples)
        if is_artifact:
            artifact_motifs.append(motif)
    return artifact_motifs


def normalize(df, normalize_factor, normalize_method_hits, normaliza_section, rank_method, fixed_min, fixed_max):
    # For factor log, otherwise it is linear:
    if normalize_factor == 'log':
        df = normalize_log(df, rank_method)
    
    min_motifs = df.min() 
    max_motifs = df.max()
    min_exp = min(min_motifs)
    max_exp = max(max_motifs)
    # min_max per motif:
    if normalize_method_hits == 'min_max' and normaliza_section == 'per_motif':
        normalized_df = normalize_max_min(df, max_motifs, min_motifs)
    # min_max per exp:
    if normalize_method_hits =='min_max' and normaliza_section == 'per_exp':
        normalized_df = normalize_max_min(df, max_exp, min_exp)
    # max per motif:
    if normalize_method_hits =='max' and normaliza_section == 'per_motif':
        normalized_df = normalize_max(df, max_motifs)
    # max per motif:
    if normalize_method_hits =='max' and normaliza_section == 'per_exp':
        normalized_df = normalize_max(df, max_exp)
    # fixed_min_max:
    if normalize_method_hits =='fixed_min_max':
        if fixed_max is None:
            fixed_max = max_exp
        if fixed_min is None:
            fixed_min = 0
        df = min_max_fixed(df, fixed_max, fixed_min)
        normalized_df = normalize_max_min(df, fixed_max, fixed_min)
    return normalized_df


def write_results(df, df_statistical, pass_motifs, output_path):
    output_file_statistical = output_path + '_statistical.csv'
    output_file = output_path + '.csv'
    df_statistical.to_csv(output_file_statistical, float_format='%.3f')
    pass_motifs.insert(0, 'label')
    pass_motifs.insert(0, 'sample_name')
    df_positive_motifs = df[pass_motifs]
    df_positive_motifs.set_index('sample_name', inplace=True)
    df_positive_motifs.to_csv(output_file)


def calculation(df, label):    
    df_mean = pd.DataFrame(df.mean()).transpose()
    df_std = pd.DataFrame(df.std()).transpose()
    df_median = pd.DataFrame(df.median()).transpose()
    df_max = pd.DataFrame(df.max()).transpose()
    df_min = pd.DataFrame(df.min()).transpose()
    cal = [f'mean_{label}', f'std_{label}', f'median_{label}', f'max_{label}', f'min_{label}']
    df_calculation = pd.concat([df_mean, df_std, df_median, df_max, df_min])
    df_calculation.insert(0, 'label', cal)
    df_calculation.set_index('label', inplace=True)
    return df_calculation


def find_positive_motifs(df, threshold_mean, threshold_std, threshold_median, min_max_difference):
    positive_motifs = []
    motifs_value = []
    for motif_name in df.columns:
        if (threshold_mean is None or df.loc['mean_BC', motif_name] - df.loc['mean_other', motif_name] > threshold_mean) \
            and (threshold_std is None or df.loc['std_BC', motif_name] - df.loc['std_other', motif_name] > threshold_std) \
            and (threshold_median is None or df.loc['median_BC', motif_name] - df.loc['median_other', motif_name] > threshold_median) \
            and (not min_max_difference or df.loc['min_BC', motif_name] > df.loc['max_other', motif_name]):
            positive_motifs.append(motif_name)
            motifs_value.append('positive')
        else:
            motifs_value.append('negative')

    df.loc['values'] = motifs_value
    threshold = [threshold_mean, threshold_std, threshold_median, None, None,
                 threshold_mean, threshold_std, threshold_median, None, None, None]
    df.insert(0, 'threshold', threshold)
    return df,positive_motifs


def statistical_calculation(df, output_path, done_path, invalid_mix, threshold_mean, threshold_std, threshold_median,
                            min_max_difference, rank_method, normalize_factor, normalize_method_hits, normalize_section,
                            fixed_min, fixed_max, argv):
    source_df = df.copy()
    labels = list(set(df['label']))
    biological_condition = labels[1] if labels[0]=='other' else labels[0]
    if invalid_mix:
        artifuct_motifs = is_artifact(df, biological_condition, invalid_mix)
        df.drop(artifuct_motifs, axis=1)
    source_df = df.copy()
    df = df.drop('sample_name', axis=1)
    df.set_index('label', inplace=True)
    if rank_method == 'pval':
        df = 1-df
    if normalize_factor == 'log' or  rank_method =='hits':
        df = normalize(df, normalize_factor, normalize_method_hits, normalize_section, rank_method,  fixed_min, fixed_max) 
    df_BC = df.loc[biological_condition]
    df_other = df.loc['other']
    # Calculate the statistical functions - mean, std, median 
    df_BC_statistical = calculation(df_BC, 'BC')
    df_other_statistical = calculation(df_other, 'other')
    # Concat two dataframe to one
    df_statistical = pd.concat([df_BC_statistical, df_other_statistical])
    # Left only the positive motifs
    df_statistical, positive_motifs = find_positive_motifs(df_statistical, threshold_mean, threshold_std, threshold_median, min_max_difference)
    # Write the results 
    write_results(source_df, df_statistical, positive_motifs, output_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')

if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A csv file with data matrix to model')
    parser.add_argument('output_path', type=str, help='Path to base name file for output the results')
    parser.add_argument('done_file_path', type=str, help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--invalid_mix',type=str, default=None, help='A argument to know if there is compare to naive')
    parser.add_argument('--threshold_mean', default=None,
                        type=lambda x: float(x) if 0 < float(x) < 1 
                        else parser.error(f'The threshold of the mean diffrence should be between 0 to 1'),
                        help='If the diffrenece between the mean of BC to the mean of other is bigger than the motif is seperate')
    parser.add_argument('--threshold_std', default=None, 
                        type=lambda x: float(x) if 0 < float(x) < 1 
                        else parser.error(f'The threshold of the std diffrence should be between 0 to 1'),
                        help='If the diffrenece between the std of BC to the std of other is bigger than the motif is seperate')
    parser.add_argument('--threshold_median', default=None,
                        type=lambda x: float(x) if 0 < float(x) < 1 
                        else parser.error(f'The threshold of the median diffrence should be between 0 to 1'),
                        help='If the diffrenece between the median of BC to the median of other is bigger than the motif is seperate')
    parser.add_argument('--min_max_difference', action='store_true', help='motifs is positive if the minmal val of bc is bigger than the maximal value of other')
    parser.add_argument('--rank_method', choices=['pval', 'tfidf', 'shuffles', 'hits'], default='hits', help='Motifs ranking method')
    parser.add_argument('--normalize_factor', choices=['linear', 'log'], default='linear', help='Type of factor on number for highlight them')
    parser.add_argument('--normalize_method_hits', choices=['min_max', 'max', 'fixed_min_max'], default='min_max', 
                        help='Type of method to do the normaliztion on hits data, change the values to be between 0 to 1')
    parser.add_argument('--normalize_section', choices=['per_motif','per_exp'],  default='per_motif', help='Calculate the min and max per motifs or over all the exp data')
    parser.add_argument('--fixed_min', type=int, default=None, help='In case of fixed_min_max for normalize_method_hits set the minimum value')
    parser.add_argument('--fixed_max', type=int, default=None, help='In case of fixed_min_max for normalize_method_hits set the maximum value')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    df = pd.read_csv(args.data_path)
    statistical_calculation(df, args.output_path, args.done_file_path, args.invalid_mix,
                            args.threshold_mean, args.threshold_std, args.threshold_median, args.min_max_difference,
                            args.rank_method, args.normalize_factor, args.normalize_method_hits, args.normalize_section,
                            args.fixed_min, args.fixed_max, sys.argv)
                            