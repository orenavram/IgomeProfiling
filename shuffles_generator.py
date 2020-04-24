from random import shuffle, choice, seed, random

RANDOM_SEED = 1234
DISTANCE_FACTOR = 1
STAY_AT_SAME_PLACE_PENALTY = -1
MIN_NEIGHBORS = 1
NEIGHBORS_PENALTY = -1
NEIGHBORS_POSITIVE = 0
SET_POSITION_MOVE_PROB = 1
MAX_STEPS_FACTOR = 10
CUTOFF_SCORE_FACTOR = 1
CUTOFF_SCORE_METHOD = 'subtract' # subtract or multiply

def create_board(length):
    unset = []
    positions = {}
    for i in range(length):
        unset.append(i)
        positions[i] = i
    return unset, positions


def calculate_score(positions):
    score = 0
    for i in positions:
        if positions[i] == i:
            score += STAY_AT_SAME_PLACE_PENALTY
        else:
            score += abs(i - positions[i]) * DISTANCE_FACTOR
        if (i+1) in positions:
            if positions[i+1] > positions[i] and \
                positions[i+1] - positions[i] <= MIN_NEIGHBORS:
                score += NEIGHBORS_PENALTY
            else:
                score += NEIGHBORS_POSITIVE
    return score


def is_set_move(unset):
    if len(unset) == 0:
        return False
    if SET_POSITION_MOVE_PROB == 1:
        return True
    return random() < SET_POSITION_MOVE_PROB


def listToDict(list):
    result = {}
    for i in range(len(list)):
        result[i] = list[i]
    return result


def calculate_cutoff_score(length):
    score = 0
    pos = length - 1
    for _ in range(length - 1,  length // 2 - 1, -1):
        score += pos
        pos -= 2
    score *= 2
    if CUTOFF_SCORE_METHOD == 'multiply':
        cutoff = score * CUTOFF_SCORE_FACTOR
    else:
        cutoff = score - CUTOFF_SCORE_FACTOR
    return score, cutoff


def swap(positions, index, new_index):
    a = positions[index]
    b = positions[new_index]
    positions[index] = b
    positions[new_index] = a


def get_best_swaps(index, positions):
    scores = {}
    length = len(positions)
    for i in range(length):
        new_positions = positions.copy()
        swap(new_positions, index, i)
        score = calculate_score(new_positions)
        try:
            scores[score].append(new_positions)
        except:
            scores[score] = [new_positions]
    max_score = max(scores.keys())
    return scores[max_score]


def play(unset, positions, cutoff):
    length = len(unset)
    indexes = [i for i in range(length)]
    for _ in range(length * MAX_STEPS_FACTOR):
        if is_set_move(unset):
            index = choice(unset)
            unset.remove(index)
        else:
            index = choice(indexes)
        best_positions = get_best_swaps(index, positions)
        positions = choice(best_positions)
        score = calculate_score(positions)
        if score >= cutoff:
            break
        

    score = calculate_score(positions)
    return positions, score


def max_shuffle(length, iterations=1, file=None, output_hpp=False):
    unset, positions = create_board(length)
    # print(calculate_score(positions))
    seed(RANDOM_SEED)
    max_score, cutoff = calculate_cutoff_score(length)
    
    if file:
        def writeLine(*args):
            if output_hpp:
                if len(args) == 1:
                    file.write(f'{args[0]}\n')    
                else:
                    file.write(f'\t\t\t{{{",".join([str(x) for x in args[0]])}}},\n')
            else:
                file.write(' '.join([str(arg) for arg in args]))
                file.write('\n')
        print_func = writeLine
    else:
        print_func = print

    results = {}
    for _ in range(iterations):
        iteration_positions, iteration_score = play(list(unset), list(positions), cutoff)
        results[tuple(iteration_positions)] = iteration_score
    if output_hpp:
        print_func(f'\t{{\t{length}, {{')
    else:
        print_func(f'length={length}, max_score={max_score}, cutoff={cutoff}')
    for result in results:
        print_func(result, results[result])
    if output_hpp:
        print_func(f'\t\t}}\n\t}},')


if __name__ == '__main__':
    min_length = 1
    max_length = 20
    iterations = 100
    output_hpp = True
    output_ext = 'hpp' if output_hpp else 'txt'
    
    with open(f'controlled_shuffles.{output_ext}', 'w') as f:
        if output_hpp:
            f.write('#include <map>\n')
            f.write('#include <vector>\n')
            f.write('using namespace std;\n')
            f.write('\n')
            f.write('typedef vector<int> ShufflePattern;\n')
            f.write('typedef vector<ShufflePattern> ShufflePatterns;\n')
            f.write('typedef map<int, ShufflePatterns> ShufflesMap;\n')
            f.write('ShufflesMap _shuffle_sequences = {\n')
        for length in range(min_length, max_length + 1):
            max_shuffle(length, iterations, f, output_hpp)
        if output_hpp:
            f.write('};\n')
