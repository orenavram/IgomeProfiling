import os
from sys import argv

path_to_features = argv[1]
if len(argv)>2:
    positive_class_types = [argv[2]]
    classes_type = 'binary'
else:
    positive_class = ''
    classes_type = 'multi'

# if positive_class=='Community_N':
#     positive_class = 'COMMU_N'
# if positive_class=='Community_P':
#     positive_class = 'COMMU_P'
#
# if positive_class=='Neutralizer':
#     positive_class = 'Neut'
# if positive_class=='Nonneutralizer':
#     positive_class = 'NN'
# if positive_class=='Healthy':
#     positive_class = 'UH'

if positive_class_types==['Pre_E5']:
    positive_class_types = ['1_1', '1_2']
if positive_class_types==['Pre_E9']:
    positive_class_types = ['1_3', '1_4']
if positive_class_types==['Pre_GST']:
    positive_class_types = ['1_5', '1_6']
if positive_class_types==['Final_E5']:
    positive_class_types = ['6_1', '6_2']
if positive_class_types==['Final_E9']:
    positive_class_types = ['6_3', '6_4']
if positive_class_types==['Final_GST']:
    positive_class_types = ['6_5', '6_6']
if positive_class_types==['Pre']:
    positive_class_types = ['1_', '1_']
if positive_class_types==['Final']:
    positive_class_types = ['6_', '6_']

print(f'Positive label is: {positive_class_types} (labels of {argv[2]}) (binary labeling)')# if positive_class else "no positive class (multi)"}')
# print(f'Positive label is: {positive_class+ " (binary)" if positive_class else "no positive class (multi)"}')
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
        if classes_type == 'binary' and not any(positive_class.lower() in line_tokens[-1].lower() for positive_class in positive_class_types):
            line_tokens.insert(-1, 'other')
        else:
            #multi classes or binary with positive class
            line_tokens.insert(-1, argv[2])

        #get rid of the long names
        line_tokens[-1] = line_tokens[-1].replace('sample_NA_all.', '').replace('.MistakeAllowed.1.AA.UNIQUE.NORMALIZATION_RpM.S100.txt', '')

        result += ','.join(line_tokens) + '\n'

prefix, suffix = os.path.splitext(path_to_features)
path_to_labeled_features = f'{prefix}_{classes_type}_labeled{suffix}'
with open(path_to_labeled_features, 'w') as f:
    f.write(result)





