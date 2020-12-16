
from auxiliaries.pipeline_auxiliaries import *

def is_same_samples(samples2bc_path,barcode2samples_path):
    #verify that file sample2biologicalcondition and file barcode2sample have the same samples.
    samples_from_s2bc=sorted(load_table_to_dict(samples2bc_path,'Samples {} belongs to more than one bc!!').keys())
    samples_from_b2s=load_table_to_dict(barcode2samples_path,'Barcode {} belongs to more than one sample!!').values()
    samples_from_b2s=sorted(samples_from_b2s)
    return samples_from_b2s==samples_from_s2bc

def is_valid_structure_file(file_path):
    #verift the structure of the file, should look like: word\tword\n
    f=open(file_path)
    lines=f.readlines()
    for line in lines:
        if line.count('\t')!=1:
            return False 
        if len(line.split('\t'))!=2:
            return False
    return True        
'''
if __name__ == '__main__':
    path_file_bc='mock_data/samplename2biologicalcondition.txt'
    path_file_barcode='mock_data/barcode2samplename.txt'
    print(is_valid_structure_file(path_file_bc))
'''    