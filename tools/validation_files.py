import sys
from auxiliaries.pipeline_auxiliaries import load_table_to_dict
import os
import re
import json
import jsonschema
from jsonschema import validate

# Describe what kind of json you expect.
sample2bc_Schema = {
    "type": "object",
    "propertyNames": {
        "pattern": "^[A-Za-z0-9_]*$"
    },
    "patternProperties": {
        "^[A-Za-z0-9_]*$" : { "type": "array", "items": {"type": "string", "pattern": "^[A-Za-z0-9_]*$"} } 
    }
}


def is_same_samples(samples2bc_path, barcode2samples_path):
    # verify that file sample2biologicalcondition and file barcode2sample have the same samples.
    samples_from_sample2bc = sorted(load_table_to_dict(samples2bc_path,'Samples {} belongs to more than one bc!!').keys())
    samples_from_barcode2sample = sorted(load_table_to_dict(barcode2samples_path,'Barcode {} belongs to more than one sample!!').values())
    return samples_from_sample2bc == samples_from_barcode2sample

def is_valid_txt_file(file_path):
    # verify the structure of the file -  should look like: word\tword\n
    regex = '\w*\t\w*\n'
    for data in  open(file_path, "rb"):
        line = data.decode("utf-8")
        m = re.match(regex, line)
        if m == None:
            return False
    return True        

def is_valid_Json_file(jsonData):
    try:
        validate(instance=jsonData, schema=sample2bc_Schema)
    except jsonschema.exceptions.ValidationError as err:
        return False
    return True


def is_input_files_valid(samples2bc_path, barcode2samples_path):
    samples2bc_valid = False
    barcode2samples_valid = False
    # Test if the file samples2bc_path is valid. 
    if samples2bc_path:
        filename, file_extension = os.path.splitext(samples2bc_path)
        samples2bc_valid = is_valid_txt_file(samples2bc_path) if file_extension == '.txt' else is_valid_Json_file(json.load(open(samples2bc_path)))
        if  not samples2bc_valid: 
            print(f'The structure of file name {samples2bc_path} is not valid')   
    
    # Test if the file barcode2samples_valid is valid. 
    if barcode2samples_path:
        barcode2samples_valid = is_valid_txt_file(barcode2samples_path) 
        if not barcode2samples_valid: 
            print(f'The structure of file name {barcode2samples_path} is not valid')   
    
    # If there are two files - test that the two files have the same samples.
    if (samples2bc_path and barcode2samples_path):
        same_samples = is_same_samples(samples2bc_path, barcode2samples_path) if (samples2bc_valid and barcode2samples_valid) else True
        if not same_samples:
            print(f'The files: {samples2bc_path} and {barcode2samples_path} not contain the same samples!!')
            return False
        else:
            return False    
    else:
        # input only one file - return the validation of the input file.
        return samples2bc_valid or barcode2samples_valid