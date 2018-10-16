import os
from sys import argv

path_to_features = argv[1]
if len(argv)>2:
    positive_class = argv[2]
    elution = ''
    if len(argv) > 3:
        elution = argv[3]
    classes_type = 'binary'
else:
    positive_class = ''
    classes_type = 'multi'

result = ''
with open(path_to_features) as f:
    header_tokens = f.readline().rstrip().split(',')
    header_tokens.insert(-1, 'label')
    result += ','.join(header_tokens) + '\n'
    for line in f:
        line_tokens = line.rstrip().split(',')
        '''
        # old:
        # remove triplicate identifier '_a' / '_b' / '_c' from the end of the sample name
        # sample_name_without_triplicate_identifier = line_tokens[-1][:line_tokens[-1].rindex('_')]
        '''
        if classes_type == 'binary' and positive_class not in line_tokens[-1]:
            line_tokens.insert(-1, 'other')
        else:
            #multi classes or binary with positive class
            line_tokens.insert(-1, positive_class)
        result += ','.join(line_tokens) + '\n'

prefix, suffix = os.path.splitext(path_to_features)
path_to_labeled_features = f'{prefix}_{classes_type}_labeled{suffix}'
with open(path_to_labeled_features, 'w') as f:
    f.write(result)





