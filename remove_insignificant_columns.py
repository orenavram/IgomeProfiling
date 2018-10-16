import os
import pandas as pd
from sys import argv

path_to_features = argv[1]
positive_class = argv[2]

print(f'Filtering: {path_to_features}')
df = pd.read_csv(path_to_features)

# remove columns that have no significant values (<=0.05) in the positive labeled rows
df_positive_class = df[df['sample_name'].map(lambda s: positive_class in s)]  #  extract positive sample rows
mask = df_positive_class[df_positive_class<=0.05].any()  # mask of columns that have significant values
print('Mask is ready.')

df = df[mask.index[mask]]  # remove columns according to mask

prefix, suffix = os.path.splitext(path_to_features)
path_to_labeled_features = f'{prefix}_significant{suffix}'
df.to_csv(path_to_labeled_features, index=False)
print('Table was filtered.')





