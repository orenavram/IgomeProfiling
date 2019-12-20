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
        df = np.log2(df)
    else:  # p-values data
        df = -np.log2(df)
    cm = sns.clustermap(df, cmap="Blues", col_cluster=False, yticklabels=True)
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150/number_of_samples)
    cm.ax_heatmap.set_title(f"A heat-map of the significance of the top {number_of_features} discriminatory motifs")
    cm.savefig(f"{output_path}/{number_of_features}.svg", format='svg', bbox_inches="tight")


def plot_error_rate(errors, features, csv_file_path):
    plt.figure(dpi=1000)
    plt.plot(features, errors, '--o')
    plt.xscale('log')
    plt.ylim(-0.02, 1)
    plt.xlabel("features")
    plt.ylabel("error rate")
    plt.savefig(f"{os.path.splitext(csv_file_path)[0]}_error_rate.png")
    plt.close()


def train_models(csv_file_path):
    logging.info('Preparing output path...')
    csv_folder, csv_file_name = os.path.split(csv_file_path)
    csv_file_prefix = os.path.splitext(csv_file_name)[0]  # without extension
    output_path = os.path.join(csv_folder, csv_file_prefix)
    if not os.path.exists(output_path):
        logging.info('Creating output path...')
        os.makedirs(output_path)

    logging.info('Parsing data...')
    train_data, test_data = parse_data(csv_file_path)
    y = np.array(train_data['label']) # saving the (true) labels
    train_data.drop(['label'], axis=1, inplace=True)
    test_data.drop(['label'], axis=1, inplace=True)
    X = np.array(train_data)  # saving an array of the variables

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

            plot_heat_map(df, number_of_features, output_path, 'hits' in csv_file_path, number_of_samples)

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
        #logger.debug(model.predict(test_data.iloc[:, indexes[:number_of_features]]))

    plot_error_rate(errors, features, csv_file_path)


if __name__ == '__main__':
    logging.info('Parsing command line arguments...')
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str, help='A CSV file with data to model ')
    args = parser.parse_args()
    train_models(args.input_file)
