import sys
from auxiliaries.pipeline_auxiliaries import load_table_to_dict, sample2bc_schema
import os
import re
import json
import jsonschema


def is_same_samples(samples2bc_path, barcode2samples_path):
    # verify that file sample2biologicalcondition and file barcode2sample have the same samples.
    samples_from_sample2bc = sorted(load_table_to_dict(samples2bc_path,'Samples {} belongs to more than one bc!!').keys())
    samples_from_barcode2sample = sorted(load_table_to_dict(barcode2samples_path,'Barcode {} belongs to more than one sample!!').values())
    if not samples_from_sample2bc == samples_from_barcode2sample:
        print(f'The files: {samples2bc_path} and {barcode2samples_path} not contain the same samples!!')
        return False
    return True

def is_valid_txt_file(file_path):
    # verify the structure of the file -  should look like: word\tword\n
    regex = '\w*\t\w*\n'
    for data in  open(file_path, "rb"):
        line = data.decode("utf-8")
        m = re.match(regex, line)
        if m == None:
            print(f'The structure of file name {file_path} is not valid')
            return False
    return True        

def is_valid_Json_file(json_path, schema):
    json_data = json.load(open(json_path))
    try:
        jsonschema.validate(instance=json_data, schema=schema)
    except jsonschema.exceptions.ValidationError as err:
        print(f'The structure of file name {json_path} is not valid')
        return False
    return True

def is_input_files_valid(samplename2biologicalcondition_path, barcode2samplename_path):
    samples2bc_valid = False
    barcode2samples_valid = False
    # Test if the file samplename2biologicalcondition_path is valid. 
    if samplename2biologicalcondition_path:
        filename, file_extension = os.path.splitext(samplename2biologicalcondition_path)
        samples2bc_valid = is_valid_txt_file(samplename2biologicalcondition_path) if file_extension == '.txt' else is_valid_Json_file(samplename2biologicalcondition_path, sample2bc_schema)
    # Test if the file barcode2samples_valid is valid. 
    if barcode2samplename_path:
        barcode2samples_valid = is_valid_txt_file(barcode2samplename_path) 
    
    # If there are two files - test that the two files have the same samples.
    if (samplename2biologicalcondition_path and barcode2samplename_path):
        return is_same_samples(samplename2biologicalcondition_path, barcode2samplename_path) if (samples2bc_valid and barcode2samples_valid) else False
    else:
        return (samples2bc_valid or barcode2samples_valid) 