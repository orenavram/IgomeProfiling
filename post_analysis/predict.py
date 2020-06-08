import sys
import os
import argparse
import pandas as pd
import joblib
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def main(data_path, model_paths):
    data = pd.read_csv(data_path)
    model_with_perfect_classification = []
    # e.g., ..../model_fitting/17b/17b_significant_pvalues_model/best_model/Top_1_features_model.pkl
    for model_path in sorted(model_paths, key=lambda x: -int(x.split('/')[-1].split('_')[1])):
        number_of_features = int(os.path.split(model_path)[-1].split('_')[1])

        feature_names = pd.read_csv(model_path.replace('_model.pkl', '.csv')).drop(['sample_name'], axis=1).columns
        assert len(feature_names) == number_of_features, f'{len(feature_names)} != {number_of_features}'
        sample_name = data['sample_name']
        X = data.drop(['label', 'sample_name'], axis=1).reindex(feature_names, axis=1).values
        assert X.shape[0] > 1, 'No samples to predict...'
        y = data['label']
        model = joblib.load(model_path)
        model_score = model.score(X, y)
        error_rate = 1-model_score
        if error_rate > 0:
            logger.debug(model_path)
            predictions = model.predict(X)
            errors = predictions != y
            results = pd.DataFrame({'sample names': sample_name,
                                    'predictions': predictions,
                                    'true labels': y,
                                    'error': errors})
            logger.info(f'Model with *{number_of_features}* features has {sum(errors)} prediction errors')
            logger.info('\n' + results.reindex(['sample names', 'predictions', 'true labels', 'error'], axis=1).to_string() + '\n\n')
            output_path = f'{os.path.split(data_path)[0]}/Top_{number_of_features}_model_predictions.txt'
            with open(output_path, 'w') as f:
                f.write(f'Error rate: {error_rate}\n')
                f.write(f'sample name\ttrue label\tprediction\tprediction score\n')
                for i in range(len(predictions)):
                    f.write(f'{sample_name[i]}\t{y[i]}\t{predictions[i]}\t{int(errors[i])}\n')
        else:
            model_with_perfect_classification.append(str(number_of_features))

        if error_rate:
            new_file = False
            if not os.path.exists(f'{os.path.split(data_path)[0]}/error_rates.txt'):
                new_file = True
            with open(f'{os.path.split(data_path)[0]}/error_rates.txt', 'a') as f:
                if new_file:
                    f.write(f'number_of_features\terror_rate\tnum_of_errors\n')
                f.write(f'{number_of_features}\t\t\t{error_rate}\t\t{errors.sum()}\n')

        logger.debug(f'Prediction error rate is {error_rate} (for debugging)')

    with open(f'{os.path.split(data_path)[0]}/perfect_classifiers.txt', 'a') as f:
        f.write('\n'.join(model_with_perfect_classification[::-1]))


if __name__ == '__main__':

    logger.info(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    parser = argparse.ArgumentParser()
    parser.add_argument('model_path', type=lambda x: x.rstrip('/'), help='A path to a pickel file with a ML model')
    parser.add_argument('data_path', help='A path to a file with samples to predict')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    if os.path.isdir(args.model_path):
        # in case its a folder with several models
        model_paths = [f'{args.model_path}/{file_name}' for file_name in os.listdir(args.model_path) if
                       file_name.endswith('pkl')]
    else:
        model_paths = [args.model_path]

    logger.info('\n' + '#' * 70 + '\nMAKE SURE TO PREDICT WITH THE SAME FEATURES THE MODEL USES!\n' + '#' * 70)
    main(args.data_path, model_paths)

    # for bc in ['HIV', 'HCV']:  #['Herceptin', 'b12', '21c', '17b']:
    #     '''
    #     /Users/Oren/Dropbox/Projects/gershoni/IgOmeProfiling/max_number_of_cluster_members_per_bc_1000/exp12_001/analysis/model_fitting/17b/17b_significant_pvalues_model/best_model
    #     /Users/Oren/Dropbox/Projects/gershoni/IgOmeProfiling/max_number_of_cluster_members_per_bc_1000/exp12_001/test/model_fitting/17b/17b_pvalues_fixed.csv
    #     '''
    #     logger.info(f'Predicting labels for {args.data_path.replace("HIV", bc).split("IgOmeProfiling/")[-1]}...')
    #     logger.info(f'Current predictions are for {bc} class. Models with less than 100% accuracy are mentioned below.\n')
    #     current_models_path = args.model_path.replace('HIV', bc)
    #     if os.path.isdir(current_models_path):
    #         # in case its a folder with several models
    #         model_paths = [f'{current_models_path}/{file_name}' for file_name in os.listdir(current_models_path) if
    #                        file_name.endswith('pkl')]
    #     else:
    #         model_paths = [current_models_path]
    #
    #     main(args.data_path.replace('HIV', bc), model_paths)
