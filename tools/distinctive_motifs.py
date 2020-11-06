
'''
Extract ranked distinctive motifs ignoring artifacts
'''
from os import path
import pandas as pd


def get_sorted_features(feature_importance_path: str):
    with open(feature_importance_path, 'r') as f:
        features = [line.strip().split('\t') for line in f.readlines()]
    features = sorted(features, key=lambda x: x[1], reverse=True)
    return features


def is_artifact(motif: str, values: pd.DataFrame, bio_cond: str, invalid_mix: str, score: float, order: int):
    data = values[['label', 'sample_name', motif]]
    other_max = data.loc[data['label'] == 'other', motif].max()
    bc_min = data.loc[data['label'] == bio_cond, motif].min()
    artifact = bc_min < other_max
    is_perfect = bc_min > other_max
    mixed_samples = []
    is_valid_mix = True
    if artifact:
        print(f'Motif {motif} is artifact, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}')
    elif is_perfect:
        print(f'Motif {motif} is perfect, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}')
    else:
        mixed_samples = list(data.loc[(data['label'] == 'other') & (data[motif] >= bc_min), 'sample_name'])
        if invalid_mix:
            is_valid_mix = not any(invalid_mix in s for s in mixed_samples)
        if is_valid_mix:
            print(f'Motif {motif} is mixed, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}, mixes={mixed_samples}')
        else:
            print(f'Motif {motif} has invalid mix, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}, mixes={mixed_samples}')
        # TODO convert mixed samples to mixed bc inc. count
    return artifact, is_perfect, is_valid_mix, mixed_samples


def extract_distinctive_motifs(count: int, epsilon: float, feature_importance_path: str, values_path: str, hits_path: str, invalid_mix: str, min_importance_score: float):
    features = get_sorted_features(feature_importance_path)
    values = pd.read_csv(values_path)
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
    for feature in features:
        i += 1
        motif = feature[0]
        score = float(feature[1])
        if score < min_importance_score:
            break
        artifact, is_perfect, is_valid_mix, mixed_samples = is_artifact(motif, values, bio_cond, invalid_mix, score, i)
        if artifact:
            artifact_count += 1
            continue
        if not is_valid_mix:
            invalid_mix_count += 1
            continue
        if total >= count and score + epsilon < last_score:
            break
        if is_perfect:
            perfect_count += 1
        else:
            for sample in mixed_samples:
                sample_count = 0
                try:
                    sample_count = mixed_count[sample]
                except:
                    pass
                mixed_count[sample] = sample_count + 1
        distinctive_motifs.append(motif)
        last_order = i
        total += 1
        if total == count:
            last_score = score
    return distinctive_motifs, last_order, perfect_count, artifact_count, invalid_mix_count, mixed_count


if __name__ == '__main__':
    base_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting'
    count = 10
    epsilon = 0
    feature_importance_path = path.join(base_path, 'ferret_texas/ferret_texas_values_exp4_model/best_model/sorted_feature_importance.txt')
    values_path = path.join(base_path, 'ferret_texas/ferret_texas_values_exp4.csv')
    hits_path = path.join(base_path, 'ferret_texas/ferret_texas_hits_exp4.csv')
    invalid_mix = 'naive'
    min_importance_score = 1e-5
    
    motifs, last_index, perfect_count, artifact_count, invalid_mix_count, mixed_count = extract_distinctive_motifs(count, epsilon, feature_importance_path, values_path, hits_path, invalid_mix, min_importance_score)
    print(f'\nDistinctive motifs ({len(motifs)}/{last_index} tested): {motifs}')
    print(f'Perfects count: {perfect_count}')
    print(f'Mixed count: {mixed_count}')
    print(f'Filtered: Artifacts count: {artifact_count}, Invalid mixes count: {invalid_mix_count}')
