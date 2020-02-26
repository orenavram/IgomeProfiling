import sys
import os
import argparse
import seaborn as sns
import numpy as np
import pandas as pd
import joblib
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def main(data_path, model_paths):
    logger.info(f'Predicting labels for {data_path}...')
    data = pd.read_csv(data_path)
    logger.info('\n' + '#' * 100 + '\nMAKE SURE TO PREDICT WITH THE SAME FEATURES THE MODEL USES!\n' + '#' * 100)
    # e.g., ..../model_fitting/17b/17b_pvalues_model/0/Top_1_features_model.pkl
    for model_path in sorted(model_paths, key=lambda x: -int(x.split('/')[-1].split('_')[1])):
        number_of_features = int(os.path.split(model_path)[-1].split('_')[1])
        logger.info(f'Current model uses *{number_of_features}* features')

        feature_names = pd.read_csv(model_path.replace('_model.pkl', '.csv')).drop(['sample_name'], axis=1).columns
        assert len(feature_names) == number_of_features, f'{len(feature_names)} != {number_of_features}'
        sample_name = data['sample_name']
        X = data.drop(['label', 'sample_name'], axis=1).reindex(feature_names, axis=1).values
        y = data['label']
        model = joblib.load(model_path)
        model_score = model.score(X, y)
        error_rate = 1-model_score
        if error_rate > 0:
            logger.info(model_path)
            predictions = model.predict(X).tolist()
            logger.info(predictions)
            errors = predictions != y
            output_path = f'{os.path.split(model_path)[0]}/Top_{number_of_features}_model_predictions.txt'
            with open(output_path, 'w') as f:
                f.write(f'Error rate: {error_rate}\n')
                f.write(f'sample name\ttrue label\tprediction\tprediction score\n')
                for i in range(len(predictions)):
                    f.write(f'{sample_name[i]}\t{y[i]}\t{predictions[i]}\t{int(errors[i])}\n')

        if error_rate:
            with open(f'{os.path.split(model_path)[0]}/error_rates.txt', 'a') as f:
                f.write(f'number_of_features\terror_rate\tnum_of_errors\n')
                f.write(f'{number_of_features}\t\t\t{error_rate}\t\t{errors.sum()}\n')
        # if number_of_features==4:
        #     for i in range(100):
        #         if i%10==0:
        #             print(i)
        #         prev_predictions = predictions
        #         predictions = get_predictions(model, X[:, :number_of_features])
        #         if predictions != prev_predictions:
        #             print(i)
        #             print(predictions)
        #             print(prev_predictions)
        #         prev_errors = errors
        #         errors = y != predictions
        #         if (errors != prev_errors).any():
        #             print(errors)
        #             print(prev_errors)
        #         prev_error_rate = error_rate
        #         error_rate = errors.mean()
        #         if error_rate != prev_error_rate:
        #             print(error_rate)
        #             print(prev_error_rate)
        #             raise

        logger.info(f'Prediction error rate is {error_rate} (for debugging)')


if __name__ == '__main__':

    logger.info(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

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
        model_paths = [f'{args.model_path}/{file_name}' for file_name in os.listdir(args.model_path) if file_name.endswith('pkl')]
    else:
        model_paths = [args.model_path]

    main(args.data_path, model_paths)
