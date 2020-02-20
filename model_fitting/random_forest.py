import sys
import os
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.model_selection import RandomizedSearchCV
# from sklearn.metrics import plot_roc_curve

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')



def parse_data(file_path):
    data = pd.read_csv(file_path, engine='python')

    test_rows_mask = data['sample_name'].str.contains('test')
    train_rows_mask = ~test_rows_mask

    # separating train and test samples (fifth sample of each mAb)
    # a matrix of the actual feature values
    y_train = data['label'][train_rows_mask]
    y_test = data['label'][test_rows_mask]

    sample_names_train = data['sample_name'][train_rows_mask]
    sample_names_test = data['sample_name'][test_rows_mask]

    # don't need the labels and sample names anymore
    data.drop(['sample_name', 'label'], axis=1, inplace=True)
    # a matrix of the actual feature values
    X_train = data[train_rows_mask].values
    X_test = data[test_rows_mask].values

    feature_names = np.array(data.columns)

    return X_train, y_train, X_test, y_test, feature_names, sample_names_train, sample_names_test


def get_hyperparameters_grid():
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(start=100, stop=2000, num=20)]
    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']
    # Maximum number of levels in tree
    max_depth = [int(x) for x in np.linspace(10, 110, num=11)]
    max_depth.append(None)
    # Minimum number of samples required to split a node
    min_samples_split = [2, 5, 10]
    # Minimum number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4, 8]
    # Method of selecting samples for training each tree
    bootstrap = [True, False]
    # Create the random grid
    random_grid = {'n_estimators': n_estimators,
                   'max_features': max_features,
                   'max_depth': max_depth,
                   'min_samples_split': min_samples_split,
                   'min_samples_leaf': min_samples_leaf,
                   'bootstrap': bootstrap}

    # Use the random grid to search for best hyperparameters
    return random_grid


def generate_heat_map(df, number_of_features, hits, number_of_samples, output_path):
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
    raw_data_path = f"{output_path_dir}/error_rate.txt"
    with open(raw_data_path, 'w') as f:
        f.write('Features\tErrors\n')
        f.write('\n'.join(f'{feature}\t{error}' for feature, error in zip(features[::-1], errors[::-1])))
    plt.figure(dpi=1000)
    plt.plot(features, errors, '--o')
    plt.xscale('log')
    plt.ylim(-0.02, 1)
    plt.xlabel("features")
    plt.ylabel("error rate")
    plot_path = raw_data_path.replace('.txt', '.png')
    plt.savefig(plot_path)
    plt.close()


def train(X, y, hyperparameters_dict, feature_names, sample_names, hits_data, output_path):
    logger.debug('\n'+'#'*100 + f'\nTrue labels:\n{y.tolist()}\n' + '#'*100 + '\n')


    logging.info(f'Current hyper-parameters configuration is:\n{hyperparameters_dict}')
    logging.info('Training...')
    # Fit the best model configuration on the WHOLE dataset
    rf = RandomForestClassifier(**hyperparameters_dict)
    model = rf.fit(X, y)
    importance = model.feature_importances_

    # the permutation needed to get the feature importances in a decreasing order
    decreasing_feature_importance = np.argsort(importance)[::-1]
    assert (sorted(importance, reverse=True) == importance[decreasing_feature_importance]).all()

    # sort data by feature importance permutation
    importance = importance[decreasing_feature_importance]
    feature_names = feature_names[decreasing_feature_importance]
    X = X[:, decreasing_feature_importance]

    with open(f'{output_path}/feature_importance.txt', 'w') as f:
        for i in range(len(importance)):
            f.write(f'{feature_names[i]}\t{importance[i]}\n')

    number_of_samples, number_of_features = X.shape
    cv_avg_error_rate = previous_cv_avg_error_rate = 1
    cv_avg_error_rates = []
    number_of_features_per_model = []
    while cv_avg_error_rate <= previous_cv_avg_error_rate and number_of_features >= 1:
        logger.info(f'Number of features is {number_of_features}')
        number_of_features_per_model.append(number_of_features)

        # save previous cv_avg_error_rate to make sure the performances do not deteriorate
        previous_cv_avg_error_rate = cv_avg_error_rate

        # predictions = model.predict(X)
        # logger.info(f'Full model error rate is {1 - (predictions == y).mean()}')
        # logger.info(f'Current model\'s predictions\n{predictions.tolist()}')

        # compute current model accuracy for each fold of the cross validation
        cv_score = cross_val_score(model, X, y, cv=StratifiedKFold(n_splits=4), n_jobs=-1)

        # current model cv_avg_error_rate rate
        cv_avg_error_rate = 1 - cv_score.mean()
        cv_avg_error_rates.append(cv_avg_error_rate)
        logger.info(f'CV models average error rate is {cv_avg_error_rate}')

        # save current model features to a csv file
        df = save_model_features(X, feature_names, sample_names, number_of_features, output_path)

        generate_heat_map(df, number_of_features, hits_data, number_of_samples, output_path)

        #generate_roc_curve(X, y, rf, number_of_features, output_path)

        # save the model itself (serialized) for future use
        joblib.dump(model, os.path.join(output_path, f'Top_{number_of_features}_features_model.pkl'))

        # Sanity check for debugging: predicting the test data
        predictions = model.predict(X)
        logger.info(f'Full model error rate is {1 - (predictions == y).mean()} (for debugging)')
        logger.debug(f'Current model\'s predictions\n{predictions.tolist()}')

        # update number of features
        number_of_features //= 2

        # extract only the (new) half most important features
        feature_names = feature_names[:number_of_features]
        X = X[:, :number_of_features]

        if number_of_features > 0:
            # re-evaluate
            # rf = RandomForestClassifier(n_estimators=number_of_trees)
            model = rf.fit(X, y)

    return cv_avg_error_rates, number_of_features_per_model


def save_model_features(X, feature_names, sample_names, number_of_features, output_path):
    df = pd.DataFrame()
    for i in range(len(feature_names)):
        df[feature_names[i]] = X[:, i]  # train_data.iloc[:, :number_of_features]
    df.set_index(sample_names, inplace=True)
    df.to_csv(f"{output_path}/Top_{number_of_features}_features.csv")
    return df


def generate_roc_curve(X, y, classifier, number_of_features, output_path):
    ax = plt.gca()
    plot_roc_curve(classifier, X, y, ax=ax)
    plt.savefig(f"{output_path}/roc_curve_{number_of_features}.png", format='png', bbox_inches="tight")
    plt.close()


def train_models(csv_file_path, num_of_iterations, done_path, cv_num_of_splits, argv):
    logging.info('Preparing output path...')
    csv_folder, csv_file_name = os.path.split(csv_file_path)
    csv_file_prefix = os.path.splitext(csv_file_name)[0]  # without extension
    output_path = os.path.join(csv_folder, f'{csv_file_prefix}_model')
    feature_selection_summary_path = f'{output_path}/feature_selection_summary.txt'

    logging.info('Parsing data...')
    is_hits_data = 'hits' in os.path.split(csv_file_path)[-1]  # Does the file name contain "hits"?

    X_train, y_train, X_test, y_test, feature_names, sample_names_train, sample_names_test = parse_data(csv_file_path)

    logger.info('\n'+'#'*100 + f'\nTrue labels:\n{y_train.tolist()}\n' + '#'*100 + '\n')

    logging.info('Sampling hyperparameters...')
    hyperparameters_grid = get_hyperparameters_grid()

    if cv_num_of_splits < 2:
        logging.info('Number of CV folds is less than 2. Using Leave One Out approach!')
    else:
        logging.info(f'Number of CV folds is {cv_num_of_splits}.')

    # search across <n_iter> different combinations, and use all available cores
    rf_random = RandomizedSearchCV(RandomForestClassifier(), n_jobs=-1,
                                   param_distributions=hyperparameters_grid,
                                   n_iter=num_of_iterations, verbose=2, random_state=42,
                                   cv=StratifiedKFold(n_splits=cv_num_of_splits)) # same classes proportion in each fold

    # Harnessing RandomizedSearchCV to get uniformly sampled configurations
    # not that efficient but, hey, no bugs :)
    configurations = rf_random.fit(X_train, y_train)
    logging.info(f'According to RandomizedSearchCV, BEST hyper-parameter configuration is:\n{configurations.best_params_}')

    # each element in this iterated collection is a dictionary with a different hyper-parameter configurations that
    # were sampled from <hyperparameters_grid>
    for i, sampled_configuration in enumerate(configurations.cv_results_['params']):
        output_path_i = os.path.join(output_path, str(i).zfill(len(str(num_of_iterations))))
        if not os.path.exists(output_path_i):
            logging.info('Creating output path...')
            os.makedirs(output_path_i)

        save_configuration_to_txt_file(sampled_configuration, output_path_i)

        errors, features = train(X_train, y_train, sampled_configuration, feature_names, sample_names_train, is_hits_data, output_path_i)

        plot_error_rate(errors, features, output_path_i)

        with open(feature_selection_summary_path, 'w' if i == 0 else 'a') as f:  # override only on first iteration
            f.write(f'{i}\t{features[-1]}\n')

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


def save_configuration_to_txt_file(sampled_configuration, output_path_i):
    with open(f'{output_path_i}/hyperparameters_configuration.txt', 'w') as f:
        for key in sampled_configuration:
            f.write(f'{key}={sampled_configuration[key]}\n')


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A csv file with data matrix to model ')
    parser.add_argument('num_of_iterations', type=int, help='How many random configurations of hyperparameters should be sampled?')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--cv_num_of_splits', default=4, help='How folds should be in the cross validation process? (use 0 for leave one out)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    train_models(args.data_path, args.num_of_iterations, args.done_file_path, args.cv_num_of_splits, argv=sys.argv)

