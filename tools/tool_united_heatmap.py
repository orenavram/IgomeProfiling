import sys
import pandas as pd
import os
import numpy as np
import logging
import seaborn as sns
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import *

def generate_color(mount_bc):
    color = []
    for i in range(0, len(mount_bc)+1):
        color.append(tuple(np.random.choice(range(0, 2), size=3)))
    return color

def united_csv(input_path, output_path, samplename2biologicalcondition_path, verbose):
    samplename2biologicalcondition = load_table_to_dict(samplename2biologicalcondition_path, 'Barcode {} belongs to more than one sample_name!!')
    biological_conditions = sorted(set(samplename2biologicalcondition.values()))
    dic = {}
    list_hits = []
    list_values = []
    sample_name = []
    flag_firs_time = True
    for bc in biological_conditions:
        bc_dir_path = os.path.join(input_path, bc)
        bc_hits = os.path.join(bc_dir_path,f'{bc}_hits.csv')
        bc_values = os.path.join(bc_dir_path,f'{bc}_values.csv')
        df_hits = pd.read_csv(bc_hits)
        df_value = pd.read_csv(bc_values)
        dic[bc] = df_hits.columns
        if flag_firs_time:
            sample_name = df_hits['sample_name']
            flag_firs_time = False
        #remove the sample name and label.
        df_hits = df_hits.drop(['sample_name', 'label'],axis=1)
        df_value = df_value.drop(['sample_name', 'label'],axis=1)
        dic[bc] = df_hits.columns
        list_hits.append(df_hits)
        list_values.append(df_value)
    #merge all the df to one df. 
    hits_all=pd.concat(list_hits, axis=1) 
    values_all=pd.concat(list_values,axis=1)
    #add columns sample name and label
    loc=0
    column_sample_name = 'sample_name'
    #column_label = 'label'
    #label = []
    #for i in sample_name:
    #    label.append(samplename2biologicalcondition[i])
    #hits_all.insert(loc, column_label, label)
    #values_all.insert(loc, column_label, label)
    hits_all.insert(loc, column_sample_name, sample_name)
    values_all.insert(loc, column_sample_name, sample_name)
    #save as csv file
    #output_path_hits=os.path.join(output_path,'hits.csv')
    #output_path_values=os.path.join(output_path,'values.csv')
    #hits_all.to_csv(output_path_hits, index=None)
    #values_all.to_csv(output_path_values, index=None)

    return hits_all, values_all



def generate_heat_map(df, hits_data, number_of_samples, output_path):
    df = df.set_index(df.columns[0])
    train_data = np.log2(df+1) if hits_data else df 
    cm = sns.clustermap(train_data, cmap="Blues", col_cluster=False, yticklabels=True)
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150/number_of_samples)
    cm.ax_heatmap.set_title(f"A heat-map of the all biological condition discriminatory motifs")
    cm.savefig(f"{output_path}.svg", format='svg', bbox_inches="tight")
    plt.close()

if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A path for all the csv files with data matrix to model ')
    parser.add_argument('output_path', type=str, help='A path for folder to put the output')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='Path for a file samplename2biologicalcondition')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    hits,values = united_csv(args.data_path, args.output_path, args.samplename2biologicalcondition_path, args.verbose)
    
    output_path_hits = os.path.join(args.output_path,'hits_all_bc')
    output_path_values = os.path.join(args.output_path,'values_all_bc')

    generate_heat_map(hits, True, len(hits), output_path_hits)
    generate_heat_map(values, False, len(values), output_path_values)