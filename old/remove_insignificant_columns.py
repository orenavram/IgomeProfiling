import os
import pandas as pd
from sys import argv

path_to_features = argv[1]
positive_class_types = [argv[2]]

# if positive_class_types==['Community_N']:
#     positive_class_types = ['COMMU_N']
# if positive_class_types==['Community_P']:
#     positive_class_types = ['COMMU_P']

# if positive_class=='Neutralizer':
#     positive_class = 'Neut'
# if positive_class=='Nonneutralizer':
#     positive_class = 'NN'
# if positive_class=='Healthy':
#     positive_class = 'UH'

if positive_class_types==['Pre_E5']:
    positive_class_types = ['1_1', '1_2']
if positive_class_types==['Pre_E9']:
    positive_class_types = ['1_3', '1_4']
if positive_class_types==['Pre_GST']:
    positive_class_types = ['1_5', '1_6']
if positive_class_types==['Final_E5']:
    positive_class_types = ['6_1', '6_2']
if positive_class_types==['Final_E9']:
    positive_class_types = ['6_3', '6_4']
if positive_class_types==['Final_GST']:
    positive_class_types = ['6_5', '6_6']
if positive_class_types==['Pre']:
    positive_class_types = ['1_', '1_']
if positive_class_types==['Final']:
    positive_class_types = ['6_', '6_']

print(f'Filtering: {path_to_features}.\nPositive class types is: {positive_class_types} (labels of {argv[2]})')
df = pd.read_csv(path_to_features)

# remove columns that have no significant values (<=0.05) in the positive labeled rows
df_positive_class = df[df['sample_name'].map(lambda s: any(positive_class.lower() in s.lower() for positive_class in positive_class_types))]  #  extract positive sample rows
# print("df_positive_class:")
# print(df_positive_class.to_string())
mask = df_positive_class[df_positive_class<=0.05].any()  # mask of columns that have significant values
# print("mask")
# print(mask)
print(f'Mask is ready. {mask.__neg__().sum()} out of {mask.shape[0]} columns are going to be removed.')

df = df[mask.index[mask]]  # remove columns according to mask

if df.shape[1] < 3:
    raise ValueError(f"Only {len(df.columns)} were left after filtration!! (There's a bug somewhere...)")

prefix, suffix = os.path.splitext(path_to_features)
path_to_labeled_features = f'{prefix}_significant{suffix}'
df.to_csv(path_to_labeled_features, index=False)
print('Table was filtered.')





