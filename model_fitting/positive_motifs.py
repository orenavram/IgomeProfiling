import logging
import os
import sys
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
    return df / max_value


def normalize_max_min(df, max_value, min_value):
    return (df - min_value) / (max_value - min_value)


def normalize_log(df):
    return np.log2(df + 1)


def min_max_fixed(df, max_value, min_value):
    df[df > max_value] = max_value
    df[df < min_value] = min_value
    return df


def normalize(df, normalize_factor_hits, normalize_method_hits, normaliza_section):
    # For factor log, otherwise it is linear:
    if normalize_factor_hits == 'log':
        df = normalize_log(df)
    
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
        fixed_min = 0
        fixed_max = max_exp
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
    # BC: 0 - mean, 1 - std, 2 - median, 3 - max ,4 -min 
    # other: 5 - mean, 6 - std, 7 - median, 8 - max, 9 - min
    for (motif_name, motif_data) in df.iteritems():
        if ((threshold_mean != 0.0 and  motif_data[0] - motif_data[5] > threshold_mean) or threshold_mean == 0.0) \
            and ((threshold_std != 0.0 and  motif_data[1] - motif_data[6] > threshold_std) or threshold_std == 0.0) \
            and ((threshold_median != 0.0 and  motif_data[2] - motif_data[7] > threshold_median) or threshold_median ==0.0) \
            and (min_max_difference and motif_data[4] > motif_data[8] or not min_max_difference):
            positive_motifs.append(motif_name)
    return positive_motifs


def statistical_calculation(df, biological_condition, output_path, done_path, threshold_mean, threshold_std, threshold_median,
                            min_max_difference, rank_method, normalize_factor_hits, normalize_method_hits, normaliza_section, argv):
    source_df = df.copy()
    df = df.drop('sample_name', axis=1)
    df.set_index('label', inplace=True)
    if rank_method == 'hits':
        df = normalize(df, normalize_factor_hits, normalize_method_hits, normaliza_section) 
    if rank_method == 'pval':
        df = 1-df    
    df_BC = df.loc[biological_condition]
    df_other = df.loc['other']
    # Calculate the statistical functions - mean, std, median 
    df_BC_statistical = calculation(df_BC, 'BC')
    df_other_statistical = calculation(df_other, 'other')
    # Concat two dataframe to one
    df_statistical = pd.concat([df_BC_statistical, df_other_statistical])
    # Left only the positive motifs
    positive_motifs = find_positive_motifs(df_statistical, threshold_mean, threshold_std, threshold_median, min_max_difference)
    # Write the results 
    write_results(source_df, df_statistical, positive_motifs, output_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')

if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A csv file with data matrix to model')
    parser.add_argument('biological_condition', type=str, help='Positive class\' label. All samples of another biological condition will be labeled as "other"')
    parser.add_argument('output_path', type=str, help='Path to base name file for output the results')
    parser.add_argument('done_file_path', type=str, help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--threshold_mean', type=float, default=0.0, help='If the diffrenece between the mean of BC to the mean of other is smaller than the threshold the motif is positive')
    parser.add_argument('--threshold_std', type=float, default=0.0, help='If the diffrenece between the std of BC to the std of other is smaller than the threshold the motif is positive')
    parser.add_argument('--threshold_median', type=float, default=0.0, help='If the diffrenece between the median of BC to the median of other is smaller than the threshold the motif is positive')
    parser.add_argument('--min_max_difference', action='store_true', help='motifs is positive if the minmal val of bc is bigger than the maximal value of other')
    parser.add_argument('--rank_method', choices=['pval', 'tfidf', 'shuffles', 'hits'], default='hits', help='Motifs ranking method')
    parser.add_argument('--normalize_factor_hits', choices=['linear', 'log'], default='linear', help='Type of factor on number for highlight them')
    parser.add_argument('--normalize_method_hits', choices=['min_max', 'max', 'fixed_min_max'], default='min_max', 
                        help='Type of method to do the normaliztion on hits data, change the values to be between 0 to 1')
    parser.add_argument('--normaliza_section', choices=['per_motif','per_exp'],  default='per_motif', help='Calculate the min and max per motifs or over all the exp data')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    df = pd.read_csv(args.data_path)
    statistical_calculation(df, args.biological_condition, args.output_path, args.done_file_path,
                            args.threshold_mean, args.threshold_std, args.threshold_median, args.min_max_difference,
                            args.rank_method, args.normalize_factor_hits, args.normalize_method_hits, args.normaliza_section, sys.argv)
