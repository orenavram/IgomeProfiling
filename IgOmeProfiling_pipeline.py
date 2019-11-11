import os
import subprocess
import time
import datetime
import logging
logger = logging.getLogger('main')

from global_params import *

#/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py or.cmds ./ -q "pupkolab -l nodes=compute-0-298"


pssm_score_peptide = "/groups/pupko/orenavr2/gershoni/src/PSSM_score_Peptide_Jan2018.NonVerbose/PSSM_score_Peptide.verbose"
python = "python/python-anaconda3.6.5"
gcc = "gcc/gcc-8.2.0"
#que = "pupkoweb"
qsub = f"python {src_dir}/q_submitter_power.py"

input_path = args.input_path.rstrip('/')
output_path = args.output_path.rstrip('/')
first_phase_path = f'{output_path}/first_phase_output'
first_phase_logs = f'{first_phase_path}/logs'
motif_inference_path = f'{output_path}/motif_inference'
motif_inference_logs = f'{motif_inference_path}/logs'
classificaion_path = f'{output_path}/classification'
classificaion_logs = f'{classificaion_path}/logs'

barcode2samplename_path = f'{input_path}/metadata/barcode2samplename.txt'
samplename2biologicalcondition_path = f'{input_path}/metadata/samplename2biologicalcondition.txt'

os.makedirs(output_path, exist_ok=True)

os.makedirs(first_phase_path, exist_ok=True)
os.makedirs(f'{first_phase_path}/logs', exist_ok=True)

os.makedirs(motif_inference_path, exist_ok=True)
os.makedirs(f'{motif_inference_path}/logs', exist_ok=True)

os.makedirs(classificaion_path, exist_ok=True)
os.makedirs(classificaion_logs, exist_ok=True)
os.makedirs(f'{classificaion_logs}/done', exist_ok=True)
os.makedirs(f'{classificaion_path}/pvalues', exist_ok=True)
os.makedirs(f'{classificaion_path}/hits', exist_ok=True)

configuration_file_path = f'{output_path}/configuration.txt'
logger.info(f'Writing configuration file to: {configuration_file_path}')
#TODO: write configuration file properly
with open(configuration_file_path, 'w') as f:
    with open("/groups/pupko/bialik/py/new/global_params.py", 'r') as f2:
        for l in f2:
            f.write(l)

comma_separated_barcodes = ""

# logger.info(f'{datetime.datetime.now()}: Loading barcode2samplename file from:\n{barcode2samplename_path}')
# with open(barcode2samplename_path) as f:
#     for line in f:
#         comma_separated_barcodes += line.split('\t')[0] + ','
# comma_separated_barcodes = comma_separated_barcodes[:-1]  # remove last redundant comma


logger.info(f'{datetime.datetime.now()}: Loading samplename2biologicalcondition file from:\n{samplename2biologicalcondition_path}')
biological_conditions = set()
with open(samplename2biologicalcondition_path, 'r') as f:
    for line in f:
        bc = line.rstrip().split('\t')[1]
        biological_conditions.add(bc)

logger.info(f'{datetime.datetime.now()}: Filtering the following fastq file:\n{fastq_file}')
subprocess.call(
    f"python {src_dir}/filter_{reads_type}_reads.py "  # FilterReads_and_Translate_fth1_RL_P8_LongReads.py 
    f" {fastq_file} "
    f"{output_path}/first_phase_output/ "
    f"{input_path}/metadata/barcode2samplename.txt ",
    shell=True)


#need to change it to sent the coutuniqseq_fasta.py as jobs
# print(datetime.datetime.now())
for path, dirname, filename in os.walk(folder_output + "first_phase_output"):
    for f in filename:
        if f.endswith('.faa'):
            file = path + '/' + f
            subprocess.call(f"python /groups/pupko/bialik/py/countUniqSeq_Fasta.py {file} YES RpM", shell=True)
print(datetime.datetime.now())
print("RPM has been finished")



# print(datetime.datetime.now())
logger.info(f'Inferring motifs...')
subprocess.call(
    f"python /groups/pupko/bialik/py/new/motif_inference.py {folder_input}metadata/ {folder_output}first_phase_output/ {folder_output}motif_inference"
    f" {queue_name} {max_peptide} {biggest_cluster} {biggest_cluster_sec} {1-gap_threshold}", shell=True)
# print(datetime.datetime.now())


# sequence logo
with open(f'{motif_inference_logs}/weblogo.cmds','w') as f:
    for bc in biological_conditions:
        os.makedirs(f'{motif_inference_path}/{bc}/weblogo', exist_ok=True)
        for alnfile in os.listdir(f'{motif_inference_path}/{bc}/fixed_alignments'):
            aln_file_prefix = os.path.splitext(alnfile)[0]
            f.write(f"module load {python} {gcc}!@#")
            f.write(f'{src_dir}/weblogo_generator.py {motif_inference_path}/{bc}/fixed_alignments/{alnfile} '
                    f'{motif_inference_path}/{bc}/weblogo/{aln_file_prefix}.png!@#\t{aln_file_prefix}\n')
cmd = f'{qsub} {motif_inference_logs}/weblogo.cmds {motif_inference_logs} -q {queue_name}'
subprocess.call(cmd, shell=True)


# print(datetime.datetime.now())
append_write = 'w'
jobs = 0
samples_tup = []
for path, dirname, filename in os.walk(output_path + "/first_phase_output"):
    for f in filename:
        if f.endswith("RpM.fs"):  # 'RpM.TrimmedCysLoop.fas'):
            file = path + '/' + f
            samples_tup.append((file, os.stat(file).st_size))
samples = []
for i in sorted(samples_tup, key=lambda x: x[1], reverse=True):
    samples.append(i[0])

for bc in biological_conditions:
            if "test" in bc:
                continue
            if not os.path.exists(output_path + "/classification/logs/" + bc):
                os.mkdir(output_path + "/classification/logs/" + bc)
            subprocess.call(
                f"python {src_dir}/split_meme_like_motifs_files_and_cutoffs.py {output_path}motif_inference/{bc} {split_num}", shell=True)
            for i in range(len(os.listdir(output_path + "motif_inference/" + bc + "/pssm"))):
                pssm = f"{output_path}motif_inference/{bc}/pssm/{i}.txt"
                cutoffs = f"{output_path}motif_inference/{bc}/cutoffs/{i}.txt"
                for sample in samples:                    
                    sample_name = sample.split('/')[len(sample.split('/'))-2].replace('_','-')
                    with open(output_path + "/classification/Calc_Pval_for_motifs_vs_samples.S100.cmds", append_write) as f:
                        f.write(f"module load {python} {gcc}!@#{pssm_score_peptide} -pssm {pssm}"
                                f" -pssm_cutoffs {cutoffs} -seq {sample}"
                                f" -out {output_path}classification/logs/{bc}/{bc}"
                                f"_motifs_vs_sample_{sample_name}_{i}.S100.txt -NrandPSSM 100 -CalcPSSM_Pval!@#"
                                f"touch {output_path}/classification/logs/done/{jobs}.done!@#\t{bc}_vs_{sample_name}_S100\n")
                    append_write = 'a'
                    jobs += 1

subprocess.call(
    f"{qsub} {output_path}/classification/Calc_Pval_for_motifs_vs_samples.S100.cmds {output_path}/classification/logs/ -q {queue_name}", shell=True)
for bc in biological_conditions:
    while True:
        if len(os.listdir(output_path + "/classification/logs/done")) == jobs:
            break
        else:
            print(f"{len(os.listdir(output_path + '/classification/logs/done'))} jobs done out of {jobs}")
            time.sleep(60)
# print(datetime.datetime.now())

print("PSSM_score_Peptide has been finished")

# print(datetime.datetime.now())
for bc in biological_conditions:
    subprocess.call(f"python {src_dir}/pvalues_new.py {output_path}classification/logs/{bc}"
                    f" {output_path}/classification/"
                    f" {input_path}/metadata/samplename2biologicalcondition.txt", shell=True)


# print(datetime.datetime.now())

csv_files = []
print("pvalues has been finished")
for path, dirname, files in os.walk(output_path + "classification"):
    for file in files:
        if file.endswith(".csv"):
            csv_files.append(path + '/' + file)

logger.info(f'Applying random forest...')
# print(datetime.datetime.now())
for csv in csv_files:
    subprocess.call(f"python {src_dir}/random_forest.py {csv}", shell=True)
# print(datetime.datetime.now())

biological_condition = ""
with open(input_path + "/metadata/samplename2biologicalcondition.txt", 'r') as f:
    for line in f:
        line = line.rstrip().split("\t")[1]
        if line not in biological_condition and "test" not in line:
            biological_condition += line + ','
biological_condition = biological_condition.rstrip(',')

for i in ['hits','pvalues']:
    subprocess.call(f"python {src_dir}/merge.py {output_path}classification/{i} {biological_condition}", shell=True)



motif_folder = f"{folder_output}/motif_inference/"
with open(folder_output + "motif_inference/logs/weblogo.cmds",'w') as f:
    for folder in os.listdir(motif_folder):
        if folder == "logs":
            continue
        if not os.path.exists(f"{motif_folder}/{folder}/weblogo"):
            os.mkdir(f"{motif_folder}/{folder}/weblogo")
        for file in os.listdir(f"{motif_folder}/{folder}/fixed_alignments/"):
            f.write(f"module load {python}!@#")
            f.write(f"{python_scripts}weblogo_generator.py "
                    f"{motif_folder}/{folder}/fixed_alignments/{file} "
                    f"{motif_folder}/{folder}/weblogo/{file.rstrip('aln')}pdf"
                    f"\t{file.rstrip('.aln')}\n")

cmd = f"{qsub} {motif_inference_path}/logs/weblogo.cmds {motif_inference_path}/logs -q {queue_name}"
#subprocess.call(cmd, shell=True)
