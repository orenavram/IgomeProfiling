import os
import sys
from subprocess import call, run
from time import time, sleep
import global_params
import logging

logger = logging.getLogger('main')
logging.basicConfig(level=logging.INFO)


nnk_table: {str: str} = {"CGT": "R", "CGG": "R", "AGG": "R",
                         "CTT": "L", "CTG": "L", "TTG": "L",
                         "TCT": "S", "TCG": "S", "AGT": "S",
                         "GCT": "A", "GCG": "A",
                         "GGT": "G", "GGG": "G",
                         "CCT": "P", "CCG": "P",
                         "ACT": "T", "ACG": "T",
                         "CAG": "Q",
                         "TAG": "q",
                         "GTT": "V", "GTG": "V",
                         "AAT": "N",
                         "GAT": "D",
                         "TGT": "C",
                         "GAG": "E",
                         "CAT": "H",
                         "ATT": "I",
                         "AAG": "K",
                         "ATG": "M",
                         "TTT": "F",
                         "TGG": "W",
                         "TAT": "Y"}

def verify_file_is_not_empty(file_path):
    import logging
    logger = logging.getLogger('main')
    # make sure that there are results and the file is not empty
    with open(file_path) as f:
        if len(f.read(10).strip()) == 0:
            # TODO: write error to a global error file
            msg = f'Input file is empty {file_path}'
            logger.error(msg)
            raise RuntimeError(msg)


def load_fasta_to_dict(fasta_path):
    """
    :param fasta_path: a path to a FASTA file
    :return: a dictionary that maps each header (without ">" and rstriped()) to its corresponding sequence (rstriped())
             an int that represent the number of sequences
             an int that represent the length of the alignment
    """
    header_to_sequence = {}
    with open(fasta_path) as f:
        for header in f:
            if not header.startswith('>'):
                raise TypeError('Illegal fasta file')
            # returns header without ">" !
            header_to_sequence[header[1:].rstrip()] = f.readline().rstrip()

    return header_to_sequence, len(header_to_sequence), len(header_to_sequence[header[1:].rstrip()])


def measure_time(total):
    hours = total // 3600
    minutes = (total % 3600) // 60
    seconds = total % 60
    if hours != 0:
        return f'{hours}:{minutes:02}:{seconds:02} hours'
    elif minutes != 0:
        return f'{minutes}:{seconds:02} minutes'
    else:
        return f'{seconds} seconds'


def wait_for_results(script_name, path, num_of_expected_results, error_file_path, suffix='done',
                     remove=False, time_to_wait=10, start=0):
    """
    :param script_name:
    :param path:
    :param num_of_expected_results:
    :param error_file_path:
    :param suffix:
    :param remove:
    :param time_to_wait:
    :param start:
    :return: waits until path contains num_of_expected_results $suffix files
    """
    if not start:
        start = time()
    logger.info(f'Waiting for {script_name}...\nContinues when {num_of_expected_results} results with suffix="{suffix}" will be in:\n{path}')
    if num_of_expected_results==0:
        logger.fatal(f'\n{"#"*100}\nnum_of_expected_results in {path} is {num_of_expected_results}!\nSomething went wrong in the previous step...\n{"#"*100}')
        #raise ValueError(f'\n{"#"*100}\nnum_of_expected_results is {num_of_expected_results}! Something went wrong in the previous step...\n{"#"*100}')
    total_time = 0
    i = 0
    current_num_of_results = 0
    while current_num_of_results < num_of_expected_results:
        sleep(time_to_wait)
        current_num_of_results = sum(1 for x in os.listdir(path) if x.endswith(suffix))
        jobs_left = num_of_expected_results - current_num_of_results
        total_time += time_to_wait
        i += 1
        if i % 5 == 0:  # print status every 5 cycles of $time_to_wait
            logger.info(f'\t{measure_time(total_time)} have passed since started waiting ({num_of_expected_results} - {current_num_of_results} = {jobs_left} more files are still missing)')
        assert not os.path.exists(error_file_path), f'An error occurred. For further details see: {error_file_path}'

    if remove:
        # execute(['python', '-u', '/groups/pupko/orenavr2/pipeline/RemoveDoneFiles.py', path, suffix])
        call(f'rm {path}/*{suffix}', shell=True)

    end = time()
    logger.info(f'Done waiting for:\n{script_name}\n(took {measure_time(int(end-start))}).\n')
    assert not os.path.exists(error_file_path), f'An error occurred. For further details see: {error_file_path}'


def submit_pipeline_step(script_path, params_lists, tmp_dir, job_name, queue_name, verbose, new_line_delimiter='!@#',
                         q_submitter_script_path='/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py',
                         required_modules_as_list=None, num_of_cpus=1):

    required_modules_as_str = 'python/python-anaconda3.6.5-orenavr2'
    if required_modules_as_list:
        # don't forget a space after the python module!!
        required_modules_as_str += ' ' + ' '.join(required_modules_as_list)
    cmds_as_str = f'module load {required_modules_as_str}'
    # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the rows in q_submitter
    cmds_as_str += new_line_delimiter

    for params in params_lists:
        cmds_as_str += ' '.join(['python', script_path, *params] + (['-v'] if verbose else [])) + ';'
        cmds_as_str += new_line_delimiter

    # GENERATE DONE FILE
    # write an empty string (like "touch" command)
    # cmds_as_str += ' '.join(['python', done_files_script_path, os.path.join(tmp_dir, job_name + '.done'), 'done'])+';'
    # cmds_as_str += new_line_delimiter

    cmds_as_str += '\t' + job_name + '\n'
    logger.debug(cmds_as_str)
    cmds_path = os.path.join(tmp_dir, job_name + '.cmds')
    with open(cmds_path, 'w') as f:
        f.write(cmds_as_str)

    # process_str = f'{q_submitter_script_path} {cmds_path} {tmp_dir} -q {queue_name} --cpu {num_of_cpus}'
    process = [q_submitter_script_path, cmds_path, tmp_dir, '-q', queue_name, '--cpu', str(num_of_cpus)]
    logger.info(f'Calling:\n{" ".join(process)}')
    run(process)


def fetch_cmd(script_name, parameters, verbose, error_path):
    cmd = f'python3 {script_name} ' + ' '.join(parameters + (['-v'] if verbose else []))
    logger.info(f'Executing:\n{cmd}')
    try:
        run(cmd, shell=True)
        # logger.info(f'Finished:\n{cmd}')
    except Exception as e:
        fail(error_path, e)



def load_barcode_to_sample_name(barcode2samplename_path):
    barcode_to_samplename = {}
    with open(barcode2samplename_path) as f:
        for line in f:
            if line.isspace():  # empty line
                continue
            barcode, sample_name = line.strip().split()
            if barcode in barcode_to_samplename:
                assert False, f'Barcode {barcode} belongs to more than one sample!!'  # TODO: write to a global error log file
            barcode_to_samplename[barcode] = sample_name
    return barcode_to_samplename


def fail(error_path, e):
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    with open(error_path, 'w') as f:
        f.write(f'\n{"$"*100}\n\n{fname}: {exc_type}, '
                f'at line: {exc_tb.tb_lineno}\n\ne.args: {e.args}\n\n{"$"*100}')