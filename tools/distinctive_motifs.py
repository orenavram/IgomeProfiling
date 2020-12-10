
'''
Extract ranked distinctive motifs ignoring artifacts
'''
from os import path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys 
colors_map = {
    'perfect': 'green',
    'mixed': 'blue',
    'artifact': 'red',
    'negative': 'orange'
}

def get_sorted_features(feature_importance_path: str):
    with open(feature_importance_path, 'r') as f:
        features = [line.strip().split('\t') for line in f.readlines()]
    features = sorted(features, key=lambda x: float(x[1]), reverse=True)
    return features


def is_artifact(motif: str, values: pd.DataFrame, bio_cond: str, invalid_mix: str, score: float, order: int, txt_file_out: str):
    data = values[['label', 'sample_name', motif]]
    other_max = data.loc[data['label'] == 'other', motif].max()
    bc_min = data.loc[data['label'] == bio_cond, motif].min()
    artifact = bc_min < other_max
    is_perfect = bc_min > other_max
    is_valid_mix = True
    mixed_samples = list(data.loc[(data['label'] == 'other') & (data[motif] >= bc_min), 'sample_name'])
    if artifact:
        txt_file_out.write(f'Motif {motif} is artifact, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}\n')
    elif is_perfect:
        txt_file_out.write(f'Motif {motif} is perfect, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}\n')
    else:
        if invalid_mix:
            is_valid_mix = not any(invalid_mix in s for s in mixed_samples)
        if is_valid_mix:
            txt_file_out.write(f'Motif {motif} is mixed, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}\n')
        else:
            txt_file_out.write(f'Motif {motif} is invalid mix, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}\n')
        # TODO convert mixed samples to mixed bc inc. count
    return artifact, is_perfect, is_valid_mix, mixed_samples


def generate_heatmap(base_path: str, df: pd.DataFrame, colors, title: str ,is_hit: bool):
    print('Generating heatmap...')
    if is_hit:
        sample_name=df['sample_name']
        df= df.drop(columns=['sample_name'])
        df=np.log2(df+1)
        df.insert(loc=0, column='sample_name', value=sample_name)

    df.set_index('sample_name', inplace=True)
    map_path = f'{base_path}.svg'
    number_of_samples = df.shape[0]

    map = sns.clustermap(df, cmap="Blues", col_cluster=False, yticklabels=True, col_colors=colors)
    plt.setp(map.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150 / number_of_samples)
    map.ax_heatmap.set_title(title, pad=25, fontsize=14)
    map.savefig(map_path, format='svg', bbox_inches="tight")
    plt.close()


def save_output(base_path: str, data):
    print('Saving results...')
    output_path = f'{base_path}.csv'
    df = pd.DataFrame(data, columns=['motif', 'label', 'is_artifact', 'is_perfect', 'is_valid_mix', 'mixed_samples', 'importance', 'order'])
    df.to_csv(output_path, index=False)


def extract_distinctive_motifs(count: int, epsilon: float, feature_importance_path: str, values_path: str, hits_path: str, invalid_mix: str, min_importance_score: float, output_base_path: str, heatmap_title: str):
    features = get_sorted_features(feature_importance_path)
    values = pd.read_csv(values_path)
    is_hits= True if values_path.find('hits')>0 else False
    # hits = pd.read_csv(hits_path, index_col=[1,0])
    # TODO check if motif is backed by hits (only log, no filter)
    bio_conds = list(values['label'].unique())
    bio_conds.remove('other')
    bio_cond = bio_conds[0]
    
    total = 0
    last_score = 0
    last_order = 0
    distinctive_motifs = []
    i = 0
    perfect_count = 0
    artifact_count = 0
    invalid_mix_count = 0
    mixed_count = {}
    colors = []
    output = []
    output_label = ''
    #create txt file for all the prints
    out_txt_path = f'{output_base_path}.txt'
    txt_file_out=open(out_txt_path,'w')
    for feature in features:
        i += 1
        motif = feature[0]
        score = float(feature[1])
        if score < min_importance_score:
            break
        artifact, is_perfect, is_valid_mix, mixed_samples = is_artifact(motif, values, bio_cond, invalid_mix, score, i,txt_file_out)
        if total >= count and score + epsilon < last_score:
            break
        if artifact:
            artifact_count += 1
            colors.append(colors_map['artifact'])
            output_label = 'artifact'
            output.append([motif, output_label, artifact, is_perfect, is_valid_mix, mixed_samples, score, i])
            continue
        if not is_valid_mix:
            invalid_mix_count += 1
            colors.append(colors_map['negative'])
            output_label = 'negative'
            output.append([motif, output_label, artifact, is_perfect, is_valid_mix, mixed_samples, score, i])
            continue
        if is_perfect:
            perfect_count += 1
            colors.append(colors_map['perfect'])
            output_label = 'perfect'
        else:
            colors.append(colors_map['mixed'])
            output_label = 'mixed'
            for sample in mixed_samples:
                sample_count = 0
                try:
                    sample_count = mixed_count[sample]
                except:
                    pass
                mixed_count[sample] = sample_count + 1
        distinctive_motifs.append(motif)
        output.append([motif, output_label, artifact, is_perfect, is_valid_mix, mixed_samples, score, i])
        last_order = i
        total += 1
        if total == count:
            last_score = score
    
    if output_base_path:
        #columns = ['sample_name'] + [x[0] for x in features[:i-1]]
        columns = ['sample_name'] + [x[0] for x in features[:count]]
        generate_heatmap(output_base_path, values[columns], colors, heatmap_title,is_hits)
        save_output(output_base_path, output)
    
    txt_file_out.write(f'\nDistinctive motifs ({len(distinctive_motifs)}/{last_order} tested): {distinctive_motifs}\n')
    txt_file_out.write(f'Perfects count: {perfect_count}\n')
    txt_file_out.write(f'Mixed count: {mixed_count}\n')
    txt_file_out.write(f'Filtered: Artifacts count: {artifact_count}, Invalid mixes count: {invalid_mix_count}\n')
    txt_file_out.close()


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('base_path', type=str, help='A path in which each subfolder corresponds to a results of the last level')
    # example:'/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting'
    parser.add_argument('count',type=int ,help='number of motifs at the end') #all=1000, top_10 =10
    parser.add_argument('feature_importance_path', type=str, help='A path for the feature importance path for specifice biological condition')
    # example:'Naive/Naive_values_exp9_model/best_model/sorted_feature_importance.txt'
    parser.add_argument('values_path',type=str ,help='A csv file of values')
    #exmple:'Naive/Naive_values_exp9.csv'    
    parser.add_argument('hits_path',type=str ,help='A csv file of hits')
    #exmple:'Naive/Naive_hits_exp9.csv'  
    parser.add_argument('output_base_path',type=str, help='A path for puting the results files in')
    #exmaple:'Naive/distinctive_motifs_all'
    parser.add_argument('--invalid_mix',type=str,default=None,help='A argument to know if there is compare to naive')
    parser.add_argument('--epsilon',type=int, default=0, help='range of mistake')
    parser.add_argument('--min_importance_score',type=int, default=0)
    args = parser.parse_args()

    feature_importance_path= path.join(args.base_path, args.feature_importance_path)
    values_path=path.join(args.base_path,args.values_path)
    hits_path=path.join(args.base_path,args.hits_path)
    output_base_path=path.join(args.base_path,args.output_base_path)
    heatmap_title='Exp 7/8 motifs (all) 234_N_HIV on 27/28 samples hits Distinctive Motifs (hits - top 10)' #change the title heatmap
    extract_distinctive_motifs(args.count, args.epsilon, feature_importance_path, values_path, hits_path, args.invalid_mix, args.min_importance_score, output_base_path, heatmap_title)
    
    
    