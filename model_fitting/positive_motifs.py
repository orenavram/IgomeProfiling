import logging
import os
import sys
import pandas as pd
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

def normalization(df):
    # Linear min max normalization
    normalized_df=(df-df.min())/(df.max()-df.min())
    # Linear max normalization
    # column_maxes = df.max()
    # df_max = column_maxes.max()
    # normalized_df = df / df_max
    return normalized_df

def write_results(df, df_statistical, pass_motifs, output_path):
    output_file_statistical = output_path + '_statistical.csv'
    output_file = output_path + '.csv'
    df_statistical.to_csv(output_file_statistical)
    df = df[[pass_motifs]]
    df.to_csv(output_file)

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
    for motifs in df.columns():
        print(motifs)

    return positive_motifs

def statistical_calculation(df, biological_condition, output_path, done_path,
                            threshold_mean, threshold_std, threshold_median,
                            min_max_difference, rank_method, argv):
    df = df.drop('sample_name', axis=1)
    df.set_index('label', inplace=True)
    if rank_method == 'hits':
        df = normalization(df)
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
    write_results(df,df_statistical, positive_motifs, output_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')

if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A csv file with data matrix to model ')
    parser.add_argument('biological_condition', help='Positive class\' label. All samples of another biological condition will be labeled as "other"')
    parser.add_argument('output_path', type=str, help='Path to base name file for output the results')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--threshold_mean', type=float, default=0.0, help='If the diffrenece between the mean of BC to the mean of other is smaller than the threshold the motif is positive')
    parser.add_argument('--threshold_std', type=float, default=0.0, help='If the diffrenece between the std of BC to the std of other is smaller than the threshold the motif is positive')
    parser.add_argument('--threshold_median', type=float, default=0.0, help='If the diffrenece between the median of BC to the median of other is smaller than the threshold the motif is positive')
    parser.add_argument('--min_max_difference', type=bool, help= 'motifs is positive if the minmal val of bc is bigger than the maximal value of other')
    parser.add_argument('--rank_method', choices=['pval', 'tfidf', 'shuffles', 'hits'], default='hits', help='Motifs ranking method')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    df = pd.read_csv(args.data_path)
    statistical_calculation(df, args.biological_condition, args.output_path, args.done_file_path,
                            args.threshold_mean, args.threshold_std, args.threshold_median,
                            args.min_max_difference, args.rank_method, sys.argv)
