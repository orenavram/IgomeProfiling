import json 
import sys
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, cross_validate
from sklearn.metrics import plot_roc_curve
import pandas as pd

def plot_error_rate(errors, features, cv_num_of_splits, output_path_dir):
    raw_data_path = f"{output_path_dir}/error_rate.txt"
    with open(raw_data_path, 'w') as f:
        f.write('Features\tErrors\n')
        f.write('\n'.join(f'{feature}\t{error}' for feature, error in zip(features[::-1], errors[::-1])))
    plt.figure(dpi=1000)
    plt.plot(features, errors, '--o')
    plt.xscale('log')
    plt.ylim(-0.02, 1)
    plt.xlabel('# of Features')
    plt.ylabel(f'{cv_num_of_splits}Fold CV Avg. Error Rate')
    plot_path = raw_data_path.replace('.txt', '.png')
    plt.savefig(plot_path)
    plt.close()

def train(rf, X, y, feature_names, sample_names, hits_data, use_tfidf, output_path, cv_num_of_splits):

    original_feature_names = feature_names[:]
    original_X = X[:]
    logger.debug('\n'+'#'*100 + f'\nTrue labels:\n{y.tolist()}\n' + '#'*100 + '\n')

    # Fit the best model configuration on the WHOLE dataset
    logging.info('Training...')
    model = rf.fit(X, y)
    importance = model.feature_importances_

    # the permutation needed to get the feature importances in a decreasing order
    decreasing_feature_importance = np.argsort(importance)[::-1]
    assert (sorted(importance, reverse=True) == importance[decreasing_feature_importance]).all()

    # the indexes that will be used in the next fitting process
    # At first, we start with all features. Next we will remove less important once.
    features_indexes_to_keep = range(len(feature_names))

    # write feature importance to storage
    with open(f'{output_path}/feature_importance.txt', 'w') as f:
        for i in range(len(importance)):
            f.write(f'{feature_names[i]}\t{importance[i]}\n')

    # write sorted feature importance to storage
    with open(f'{output_path}/sorted_feature_importance.txt', 'w') as f:
        for i in range(len(importance)):
            f.write(f'{feature_names[decreasing_feature_importance[i]]}\t'
                    f'{importance[decreasing_feature_importance[i]]}\n')

    number_of_samples, number_of_features = X.shape
    cv_avg_error_rate = previous_cv_avg_error_rate = 1
    cv_avg_error_rates = []
    number_of_features_per_model = []
    while cv_avg_error_rate <= previous_cv_avg_error_rate and number_of_features >= 1:

        # save previous cv_avg_error_rate to make sure the performances do not deteriorate
        previous_cv_avg_error_rate = cv_avg_error_rate

        # predictions = model.predict(X)
        # logger.info(f'Full model error rate is {1 - (predictions == y).mean()}')
        # logger.info(f'Current model\'s predictions\n{predictions.tolist()}')

        # compute current model accuracy for each fold of the cross validation
        cv_score = cross_val_score(model, X, y, cv=StratifiedKFold(cv_num_of_splits))

        # current model cv_avg_error_rate rate
        cv_avg_error_rate = 1 - cv_score.mean()
        number_of_features_per_model.append(number_of_features)
        cv_avg_error_rates.append(cv_avg_error_rate)

        logger.info(f'Number of features is {number_of_features} with avg. error rate of {cv_avg_error_rate}')
        if cv_avg_error_rate > previous_cv_avg_error_rate:
            # Stop procedure and discard current stats
            break

        # save current model (unsorted) features to a csv file
        df = save_model_features(original_X, features_indexes_to_keep, feature_names, sample_names, f'{output_path}/Top_{number_of_features}_features')

        generate_heat_map(df, number_of_features, hits_data, number_of_samples, use_tfidf, f'{output_path}/{number_of_features}')

        generate_roc_curve(X, y, rf, number_of_features, output_path)

        # save the model itself (serialized) for future use
        joblib.dump(model, os.path.join(output_path, f'Top_{number_of_features}_features_model.pkl'))

        # Sanity check for debugging: predicting the test data
        model_score = model.score(X, y)
        if 1-model_score > cv_avg_error_rate:
            predictions = model.predict(X).tolist()
            logging.error('1-model_score > cv_avg_error_rate !!!')
            logger.info(f'Full model error rate is {1-model_score}')
            logger.info(f'Current model\'s predictions\n{predictions}')
            logger.info(f'number_of_features {number_of_features}')
            logger.info(f'output_path {output_path}')

        # update number of features
        number_of_features //= 2

        # extract only the (new) half most important features
        features_indexes_to_keep = sorted(decreasing_feature_importance[:number_of_features])
        feature_names = original_feature_names[features_indexes_to_keep]
        X = original_X[:, features_indexes_to_keep]

        if number_of_features > 0:
            # re-evaluate
            model = rf.fit(X, y)

    return cv_avg_error_rates, number_of_features_per_model

def pre_train(configuration,csv_file_path,is_hits_data,output_path_i,model_number,feature_selection_summary_path,done_file_path,use_tfidf,cv_num_of_splits,argv):
    feature_selection_summary_f = open(feature_selection_summary_path, 'w')
    rf = RandomForestClassifier(**configuration)
    X_train, y_train, X_test, y_test, feature_names, sample_names_train, sample_names_test = parse_data(csv_file_path)
    errors, features = train(rf, X_train, y_train, feature_names, sample_names_train, is_hits_data, use_tfidf, output_path_i, cv_num_of_splits)
    plot_error_rate(errors, features, cv_num_of_splits, output_path_i)
    feature_selection_summary_f.write(f'{model_number}\t{features[-1]}\t{errors[-1]}\n')
    
    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')

if __name__ == '__main__':
print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('configuration_model', type=str, help='A dictionary of the configuration to this model ')
    parser.add_argument('csv_file_path',type=str,help='')
    parser.add_argument('is_hits_data',type=bool,help='if is hits data or pvalue data')
    parser.add_argument('output_path_i',type=str,help='A path for the results of this model')
    parser.add_argument('model_number',type=str,help='number of model')
    parser.add_argument('feature_selection_summary_path',type=str,help='A path for a file for write the results featurs')
    parser.add_argument('done_file_path',type=str,help='A path to a file that signals that the script finished running successfully')
    parser.add_argument('--tfidf', action='store_true', help="Are inputs from TF-IDF (avoid log(0))")
    parser.add_argument('--cv_num_of_splits',type=int,help='number of CV folds')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    #  rf = RandomForestClassifier(**{'n_estimators': 2000, 'max_features': 'auto', 'max_depth': 70, 'min_samples_split': 5, 'min_samples_leaf': 1, 'bootstrap': True, 'n_jobs': -1}
    configuration=json.load(args.configuration_model)
    pre_train(configuration,args.csv_file_path,args.is_hits_data,args.output_path_i,args.model_number,args.feature_selection_summary_path,args.done_file_path,args.use_tfidf,args.cv_num_of_splits,argv=sys.argv)