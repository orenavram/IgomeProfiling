"""
Extract identified biological conditions per motif
"""
from os import path
import pandas as pd


def load_samples_map(samples_to_bio_conds_path: str):
    with open(samples_to_bio_conds_path, 'r') as f:
        return dict([line.strip().split('\t') for line in f.readlines()])

def identify(values_path: str, samples_to_bio_conds_path: str, assumed_motifs, positive_bio_cond: str, negative_bio_cond: str):
    samples_to_bio_conds = load_samples_map(samples_to_bio_conds_path)
    values = pd.read_csv(values_path)
    values['label'] = [samples_to_bio_conds[sample] for sample in values['sample_name']]
    motifs = list(values.columns)[2:]
    
    total_positive = 0
    total_negative = 0
    total_others = 0
    total_perfect = 0
    total_per_bio_cond = {}
    assumption_in = []
    assumption_out = []
    assumption_others = []
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
        if assumed_motifs:
            is_assumed = motif in assumed_motifs
            if is_assumed and is_positive:
                assumption_in.append(motif)
            elif is_assumed:
                assumption_out.append(motif)
            elif is_positive:
                assumption_others.append(motif)
        print(f'{motif}: positive={is_positive}, negative={is_negative}, perfect={is_perfect}, identifies={identified_bio_conds}')
    
    print(f'\nTotal motifs={len(motifs)}, positives={total_positive}, negative={total_negative}, others={total_others}, perfects={total_perfect}')
    print(f'Count per biological condition: {total_per_bio_cond}')
    if assumed_motifs:
        print(f'\nAssumed count motifs={len(assumed_motifs)}, in={len(assumption_in)}, out={len(assumption_out)}, not expected={len(assumption_others)}')
        if assumption_in: print(f'Correctly assumed: {assumption_in}')
        if assumption_out: print(f'Incorrectly assumed: {assumption_out}')
        if assumption_others: print(f'Unexpected positives: {assumption_others}')


if __name__ == '__main__':
    base_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting'
    assumed_motifs = ['CSTTTTTRTC', 'CDTDGWIRMDPC', 'CPSHISLRNPC', 'CTRPTRPSLWMSPAC', 'CTNLCGTLAC', 'CCHHIFQRPPC', 'CPACSLFTC', 'CDAHSLFFPC', 'CGGTAGSFSRC', 'CLSSLFTRC'] # []
    values_path = path.join(base_path, 'ferret_texas/ferret_texas_values_exp9_complete.csv')
    hits_path = path.join(base_path, 'ferret_texas/ferret_texas_hits_exp9_complete.csv') # TODO check hits "back-up"
    samples_to_bio_conds_path = path.join(base_path, 'exp4+9_all_samples2bioconds.txt')
    positive_bio_cond = 'Texas'
    negative_bio_cond = 'Naive'

    identify(values_path, samples_to_bio_conds_path, assumed_motifs, positive_bio_cond, negative_bio_cond)