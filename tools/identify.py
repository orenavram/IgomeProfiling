"""
Extract identified biological conditions per motif
"""
from os import path
import pandas as pd


def load_samples_map(samples_to_bio_conds_path: str):
    with open(samples_to_bio_conds_path, 'r') as f:
        return dict([line.strip().split('\t') for line in f.readlines()])

def identify(values_path: str, samples_to_bio_conds_path: str, motifs, positive_bio_cond: str, negative_bio_cond: str):
    samples_to_bio_conds = load_samples_map(samples_to_bio_conds_path)
    values = pd.read_csv(values_path)
    values['label'] = [samples_to_bio_conds[sample] for sample in values['sample_name']]
    if not motifs:
        motifs = list(values.columns)[2:]
    
    total_positive = 0
    total_negative = 0
    total_others = 0
    total_perfect = 0
    total_per_bio_cond = {}
    for motif in motifs:
        is_negative = False
        is_positive = False
        is_perfect = False
        data = values[['label', motif]].groupby("label")
        scored_data = 2 * data[motif].max() + (1 - data[motif].std()) # TODO add configurable weights
        scored_data.sort_values(0, False, True)
        identified_bio_conds = []
        best_value = None
        for label, value in scored_data.iteritems():
            if not identified_bio_conds:
                best_value = value
                identified_bio_conds.append(label)
            elif value == best_value:
                identified_bio_conds.append(label)
            else:
                break
        if positive_bio_cond:
            is_positive = any(positive_bio_cond in s for s in identified_bio_conds)
        if negative_bio_cond:
            is_negative = any(negative_bio_cond in s for s in identified_bio_conds)
        if is_positive:
            total_positive += 1
            if len(identified_bio_conds) == 1: is_perfect = True
        if is_negative: total_negative += 1
        if is_perfect: total_perfect += 1
        if not is_positive and not is_negative: total_others += 1
        for bio_cond in identified_bio_conds:
            count = 0
            try:
                count = total_per_bio_cond[bio_cond]
            except:
                pass
            total_per_bio_cond[bio_cond] = count + 1
        print(f'{motif}: positive={is_positive}, negative={is_negative}, perfect={is_perfect}, identifies={identified_bio_conds}')
    
    print(f'\nTotal motifs={len(motifs)}, positives={total_positive}, negative={total_negative}, others={total_others}, perfects={total_perfect}')
    print(f'Count per biological condition: {total_per_bio_cond}')


if __name__ == '__main__':
    base_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting'
    motifs = ['CSTTTTTRTC', 'CDTDGWIRMDPC', 'CPSHISLRNPC', 'CTRPTRPSLWMSPAC', 'CTNLCGTLAC', 'CCHHIFQRPPC', 'CPACSLFTC', 'CDAHSLFFPC', 'CGGTAGSFSRC', 'CLSSLFTRC'] # []
    values_path = path.join(base_path, 'ferret_texas/ferret_texas_values_exp4.csv')
    hits_path = path.join(base_path, 'ferret_texas/ferret_texas_hits_exp4.csv') # TODO check hits "back-up"
    samples_to_bio_conds_path = path.join(base_path, 'exp4+9_all_samples2bioconds.txt')
    positive_bio_cond = 'texas'
    negative_bio_cond = 'naive'

    identify(values_path, samples_to_bio_conds_path, motifs, positive_bio_cond, negative_bio_cond)