import pandas as pd
import os
from sys import argv
from functools import reduce
random_forest_path = argv[1]
random_forest_path = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/1_ALL/random_forest/HITS'

features_tables = []
minimums = []
for biological_condition in os.listdir(random_forest_path):
    biological_condition_path = os.path.join(random_forest_path, biological_condition)
    if not os.path.isdir(biological_condition_path):
        continue

    min_num_of_features = float('inf')
    for item in os.listdir(biological_condition_path):
        if item.startswith('Top'):  # e.g., Top_15
            curr_num_of_features = int(item.split('_')[1])
            if curr_num_of_features < min_num_of_features:
                min_num_of_features = curr_num_of_features
                min_num_of_features_path = os.path.join(biological_condition_path, item)
    minimums.append(min_num_of_features)

    for file in os.listdir(min_num_of_features_path):
        if file.endswith('Features.csv'):  # e.g., mAb_c54_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_15_Features.csv
            features_tables.append(os.path.join(min_num_of_features_path, file))

features_tables = [
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_c54/Top_8/mAb_c54_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_8_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_CR8020/Top_4/mAb_CR8020_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_4_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_c108/Top_7/mAb_c108_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_7_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_c585/Top_4/mAb_c585_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_4_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_F045/Top_2/mAb_F045_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_2_Features.csv',
'/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/8_mAbs_vs_mAbs/random_forest/mAb_c3662/Top_3/mAb_c3662_against_Exp_DP_4.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv'
]
dfs = [pd.read_csv(csv, index_col=0).drop('label', 1) for csv in features_tables]  # remove first column (index column)
final_df = reduce(lambda left,right: pd.merge(left,right,on='sample_name'), dfs) #merge all dfs
final_df = final_df.reindex(list([a for a in final_df.columns if a != 'sample_name' and 'label' not in a])+['sample_name'], axis=1) #move sample_name to last column
final_df.to_csv(random_forest_path+'/merged_features_df.csv', index=False)

for item in features_tables:
    print(item.split('/')[-1])
print(minimums)
