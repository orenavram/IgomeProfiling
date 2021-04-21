import sys
from auxiliaries.pipeline_auxiliaries import *
import os

def is_same_samples(samples2bc_path,barcode2samples_path):
    #verify that file sample2biologicalcondition and file barcode2sample have the same samples.
    filename, file_extension = os.path.splitext(samples2bc_path)
    samples_from_sample2bc=[]
    if file_extension=='txt':
        samples_from_sample2bc=sorted(load_table_to_dict(samples2bc_path,'Samples {} belongs to more than one bc!!').keys())
    #else: #json file
        #TODO call to function that return dict from json file
    samples_from_barcode2sample=sorted(load_table_to_dict(barcode2samples_path,'Barcode {} belongs to more than one sample!!').values())
    return samples_from_sample2bc==samples_from_barcode2sample

def is_valid_structure_file(file_path):
    #verift the structure of the file, should look like: word\tword\n
    #word can conation: number, big or little letter,_
    with open(file_path, "rb") as f:
        b=f.read(1)
        find_tab=False
        count_tab=0
        while b:
            num=int.from_bytes(b, byteorder=sys.byteorder)
            if num==9: #is a TAB
                count_tab+=1
                find_tab=True
            elif num==10 and find_tab and count_tab==1: #10 is new line
                find_tab=False
                count_tab=0
            elif num<48 or (num>57 and num<65) or (num>90 and num<95) or (num>95 and num<97) or (num>122):
                #48-57 number, 65-90 Big letter, 97-122 small letter, 95 is _
                return False
            b=f.read(1)
    return True        

def is_validation_files(samples2bc_path, barcode2samples_path):
    samples2bc_valid=True
    filename, file_extension = os.path.splitext(samples2bc_path)
    if file_extension=='txt':
        samples2bc_valid = is_same_samples(samples2bc_path)
    barcode2samples_valid = is_same_samples(barcode2samples_path)
    if samples2bc_valid and barcode2samples_valid:
        return is_same_samples(samples2bc_path, barcode2samples_path)
    else:
        return False