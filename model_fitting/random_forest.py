import sys
import os
import argparse
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')



def parse_data(file_path):
    # reading the CSV file if it's legal
    try:
        data = pd.read_csv(file_path, engine='python')
    except Exception as e:
        exit(f'Cannot analyse data! {e}')

    # separating train and test samples (fifth sample of each mAb)
    train_data = data[~data['sample_name'].str.contains('test')]
    test_data = data[data['sample_name'].str.contains('test')]

    # set sample names as index
    train_data.set_index('sample_name', inplace=True)
    test_data.set_index('sample_name', inplace=True)

    return train_data, test_data


def plot_heat_map(df, number_of_features, output_path, hits, number_of_samples):
    #plt.figure(dpi=3000)
    # transform the data for better contrast in the visualization
    if hits:  # hits data
        df = np.log2(df+1)  # pseudo counts
        # df = df
    else:  # p-values data
        df = -np.log2(df)
    cm = sns.clustermap(df, cmap="Blues", col_cluster=False, yticklabels=True)
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150/number_of_samples)
    cm.ax_heatmap.set_title(f"A heat-map of the significance of the top {number_of_features} discriminatory motifs")
    cm.savefig(f"{output_path}/{number_of_features}.svg", format='svg', bbox_inches="tight")
    plt.close()


def plot_error_rate(errors, features, output_path_dir):
    plt.figure(dpi=1000)
    plt.plot(features, errors, '--o')
    plt.xscale('log')
    plt.ylim(-0.02, 1)
    plt.xlabel("features")
    plt.ylabel("error rate")
    plt.savefig(f"{output_path_dir}/error_rate.png")
    plt.close()


def train_models(csv_file_path, done_path, num_of_iterations, argv):
    logging.info('Parsing data...')
    train_data, test_data = parse_data(csv_file_path)
    y = np.array(train_data['label']) # saving the (true) labels
    train_data.drop(['label'], axis=1, inplace=True)
    test_data.drop(['label'], axis=1, inplace=True)
    X = np.array(train_data)  # saving an array of the variables
    is_hits_data = 'hits' in csv_file_path

    for i in range(num_of_iterations):
        logging.info('Preparing output path...')
        csv_folder, csv_file_name = os.path.split(csv_file_path)
        csv_file_prefix = os.path.splitext(csv_file_name)[0]  # without extension
        output_path = os.path.join(csv_folder, f'{csv_file_prefix}_model', str(i))
        if not os.path.exists(output_path):
            logging.info('Creating output path...')
            os.makedirs(output_path)

        errors, features = train(X, y, is_hits_data, train_data, output_path)

        plot_error_rate(errors, features, output_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


def train(X, y, hits_data, train_data, output_path):
    logging.info('Training...')
    rf = RandomForestClassifier(n_estimators=100)
    model = rf.fit(X, y)
    importances = model.feature_importances_
    indexes = list(np.argsort(-importances))  # negate to sort in decreasing order
    number_of_samples, number_of_features = X.shape
    error = previous_error = 1
    errors = []
    features = []
    while error <= previous_error and number_of_features >= 1:
        logger.info(f'Number of features is {number_of_features}')
        features.append(number_of_features)

        # save previous error to make sure the performances do not deteriorate
        previous_error = error

        # compute current model accuracy for each fold of the cross validation
        cv_score = cross_val_score(rf, X, y, cv=3, n_jobs=-1)

        # current model error rate
        error = 1 - cv_score.mean()
        errors.append(error)
        logger.info(f'Error rate is {error}')

        # save current model features to a csv file
        df = train_data.iloc[:, indexes[:number_of_features]]
        if error <= previous_error:
            df.to_csv(f"{output_path}/Top_{number_of_features}_features.csv")

            plot_heat_map(df, number_of_features, output_path, hits_data, number_of_samples)

            # save the model itself (serialized) for future use
            joblib.dump(model,
                        os.path.join(output_path, f'Top_{number_of_features}_features_model.pkl'))

        # update number of features
        number_of_features //= 2
        if number_of_features < 1:
            continue

        # extract only the (new) half most important features
        X = np.array(train_data.iloc[:, indexes[:number_of_features]])

        # re-evaluate
        rf = RandomForestClassifier(n_estimators=100)
        model = rf.fit(X, y)

        # Sanity check for debugging: predicting the test data
        # change the logging level (line 10) to logging.DEBUG to get the predictions
        # logger.debug(model.predict(test_data.iloc[:, indexes[:number_of_features]]))
    return errors, features


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A csv file with data matrix to model ')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--num_of_iterations', default=5, help='How many should the RF run?')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    train_models(args.data_path, args.done_file_path, args.num_of_iterations, argv=sys.argv)

