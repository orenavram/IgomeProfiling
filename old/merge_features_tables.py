import pandas as pd
import os
from sys import argv
from functools import reduce
#random_forest_path = argv[1]
#random_forest_path = '/groups/pupko/orenavr2/gershoni/Experiments/ExpHaiti/analyses/family_vs_active/random_forest/'
random_forest_path = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp22/analyses/hiv_vs_hcv/random_forest/'
#random_forest_path = 'D:/Dropbox/Projects/gershoni/Experiments/ExpHaiti/analyses/family_vs_active/random_forest/'
#random_forest_path = 'D:/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/10_stem_vs_nonstem/random_forest/'
order = ['HIV', 'HCV']
min_middle_max = 'max'

features_tables = []
feature_nums = {}
for biological_condition in os.listdir(random_forest_path):
    biological_condition_path = os.path.join(random_forest_path, biological_condition)
    if not os.path.isdir(biological_condition_path):
        continue

    num_of_features_to_file_path = {}
    for item in os.listdir(biological_condition_path):
        if item.startswith('Top'):  # e.g., Top_15
            curr_num_of_features = int(item.split('_')[1])
            num_of_features_to_file_path[curr_num_of_features] = os.path.join(biological_condition_path, item)

    if min_middle_max == 'min':
        num_of_features = sorted(num_of_features_to_file_path)[0] # min? middle? max?
    elif min_middle_max == 'middle':
        middle_index = min(len(num_of_features_to_file_path)//2, 2)
        num_of_features = sorted(num_of_features_to_file_path)[middle_index] # min? middle? max?
    else: # min_middle_max == 'max'
        num_of_features = sorted(num_of_features_to_file_path)[-1] # min? middle? max?
    num_of_features_path = num_of_features_to_file_path[num_of_features]
    feature_nums[biological_condition] = str(num_of_features)

    for file in os.listdir(num_of_features_path):
        if file.endswith('Features.csv'):  # e.g., mAb_c54_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_15_Features.csv
            features_tables.append(os.path.join(num_of_features_path, file))

non_ordered_features_tables = features_tables
features_tables = []
for item in order:
    for i in range(len(non_ordered_features_tables)):
        if item.upper() in non_ordered_features_tables[i].upper():
            features_tables.append(non_ordered_features_tables.pop(i))
            break


'''
features_tables = [
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_c108/Top_56/mAb_c108_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_56_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_c54/Top_60/mAb_c54_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_60_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_c585/Top_56/mAb_c585_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_56_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_CR8020/Top_30/mAb_CR8020_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_30_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_c3662/Top_50/mAb_c3662_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_50_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_F045/Top_66/mAb_F045_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_66_Features.csv'
]

'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/9_Ferrets_vs_Ferrets/random_forest/COBRA/Top_3/COBRA_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/9_Ferrets_vs_Ferrets/random_forest/texas/Top_3/texas_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/9_Ferrets_vs_Ferrets/random_forest/swiss/Top_3/swiss_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/9_Ferrets_vs_Ferrets/random_forest/sichuan/Top_3/sichuan_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/9_Ferrets_vs_Ferrets/random_forest/naive/Top_3/naive_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv'
]
'''
try:
    dfs = [pd.read_csv(csv, index_col=0).drop('label', 1) for csv in features_tables]  # remove first column (index column)
except:
    dfs = [pd.read_csv(csv) for csv in features_tables]  # remove first column (index column)


final_df = reduce(lambda left,right: pd.merge(left,right,on='sample_name'), dfs) #merge all dfs
final_df.drop_duplicates(inplace=True)
final_df = final_df.reindex(list([a for a in final_df.columns if a != 'sample_name' and 'label' not in a])+['sample_name'], axis=1) #move sample_name to last column
final_df.to_csv(f'{random_forest_path}/{min_middle_max}_merged_features_df.csv', index=False)

for item in features_tables:
    print(item.split('/')[-1])
for item in features_tables:
    print(f"'{item}',")
print(feature_nums)
print(order)
print([feature_nums[x] for x in order])


with open(f'{random_forest_path}/{min_middle_max}_merged_features_df_classes_size.csv', 'w') as f:
    f.write(','.join([feature_nums[x] for x in order])+'\n')

with open(f'{random_forest_path}/{min_middle_max}_merged_features_df_order.csv', 'w') as f:
    f.write(','.join(order)+'\n')

