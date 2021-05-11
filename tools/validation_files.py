import sys
from auxiliaries.pipeline_auxiliaries import load_table_to_dict
import os
import re
import json
import jsonschema
import datetime

def is_valid_Json_file(json_path, schema, logger):
    json_data = json.load(open(json_path))
    try:
        jsonschema.validate(instance=json_data, schema=schema)
    except jsonschema.exceptions.ValidationError as err:
        logger.error(f'The structure of file name {json_path} is not valid')
        return False
    return True


def is_same_samples(samples2bc_dict, barcode2samples_dict, samplename2biologicalcondition_path, barcode2samplename_path, logger):
    # verify that file sample2biologicalcondition and file barcode2sample have the same samples.
    samples_from_sample2bc = sorted(samples2bc_dict.keys())
    samples_from_barcode2sample = sorted(barcode2samples_dict.values())
    if not samples_from_sample2bc == samples_from_barcode2sample:
        logger.info(f'The files: {samplename2biologicalcondition_path} and {barcode2samplename_path} not contain the same samples!!')
        return False
    return True


def is_valid_data(path_file ,dict_data, logger):
    # verify the structure of the file -  should look like: word\tword\n
    regex = "^[A-Za-z0-9_]+$"
    for key, value in  dict_data.items():
        match_key = re.match(regex, key)
        match_value = re.match(regex, value)
        if match_key is None or match_value is None:
            logger.info(f'The structure of file: {path_file} is not valid')
            return False
    return True        


def is_input_files_valid(samplename2biologicalcondition_path, barcode2samplename_path, logger):
    samples2bc_valid = True
    barcode2samples_valid = True
    barcode2sample_data = {}
    sample2bc_data = {}

     # Test if the file barcode2samplename_path is valid. 
    if barcode2samplename_path:
        try:
            barcode2sample_data =  load_table_to_dict(barcode2samplename_path,'Barcode {} belongs to more than one sample!!')
        except:
            logger.error(f'{datetime.datetime.now()}: can not load the file - {barcode2samplename_path}')
            barcode2samples_valid = False
        if barcode2samples_valid:        
            barcode2samples_valid = is_valid_data(barcode2samplename_path, barcode2sample_data, logger) 
    
    # Test if the file samplename2biologicalcondition_path is valid. 
    if samplename2biologicalcondition_path:
        filename, file_extension = os.path.splitext(samplename2biologicalcondition_path)
        try:
            sample2bc_data = load_table_to_dict(samplename2biologicalcondition_path,'Samples {} belongs to more than one bc!!')
        except:
            logger.error(f'{datetime.datetime.now()}: can not load the file - {samplename2biologicalcondition_path}')
            samples2bc_valid = False
        if samples2bc_valid:    
            samples2bc_valid = is_valid_data(samplename2biologicalcondition_path, sample2bc_data, logger)
    
    if not samples2bc_valid or not barcode2samples_valid:
        return False
    if barcode2samplename_path and samplename2biologicalcondition_path:
        return is_same_samples(sample2bc_data, barcode2sample_data, samplename2biologicalcondition_path, barcode2samplename_path, logger)
    return True     
