# IgOme Profiling pipeline

## Overview
IgOme Profiling pipeline is ...

## Getting started

### OS
The pipeline runs on Linux environment only.  
For Windows the pipeline can run using [WSL and Visual Studio Code](https://code.visualstudio.com/remote-tutorials/wsl/getting-started).

Tested on Ubuntu 18.04 via Windows.

### Prerequisites
The following are required in order to run the code:  
* git, e.g.:
  ```bash
  sudo apt install git
  ```
* Python 3.6+, .e.g:
  ```bash
  sudo apt install python3 python3-pip python3-venv
  ```
* g++, e.g.:
  ```
  sudo apt install g++
  ```
* CD-HIT:
  * Install using instructions: http://weizhongli-lab.org/cd-hit/
  * Make sure cd-hit is installed globally, e.g. /usr/local/bin
  * Test by running:
    ```bash
    # From project directory
    cd-hit -h
    ```
* MAFFT:
  * Install using instructions: https://mafft.cbrc.jp/alignment/software/
  * Make sure mafft is installed globally, e.g. /usr/bin
  * Test by running:
    ```bash
    # From project directory
    mafft -h
    ```

It is advised to update the packages manager before installing, e.g.:
```bash
sudo apt update
```

## Repository
Configure git user and email (if not done yet):
```bash
git config --global user.name "FIRST_NAME LAST_NAME"
git config --global user.email "you@mail.com"
```

Clone the repository:
```bash
# From projects/workspace directory
git clone https://github.com/Webiks/IgomeProfiling.git
cd IgomeProfiling
```

### Compile
Most of the code is in Python but some of the code is in C++ and requires compilation:
* UnitePSSMs:
  ```bash
  # From project directory
  cd UnitePSSMs
  g++ *.cpp -std=c++11 -O3 -o UnitePSSMs
  ```
* PSSM_score_Peptide:
  ```bash
  # From project directory
  cd PSSM_score_Peptide
  g++ *.cpp -std=c++11 -O3 -o PSSM_score_Peptide
  ```

### Python 
Python3 should be installed in previous steps.  
Use the following to install the required packages:
* Create a virtual environment for the project:
  ```bash
  # From project directory
  python3 -m venv .venv
  source .venv/bin/activate
  ```
* Install the dependencies:
  ```bash
  pip3 install -r requirements.txt
  ```

## Running
The pipeline is modular, each step is a python file which can be run by itself by passing command line arguments.  
To run the entire pipeline the entry point is ```IgOmeProfiling_pipeline.py```.  
Run the pipeline using the following:  
```bash
python3 IgOmeProfiling_pipeline.py fastq_file_path barcodes_to_samples_path samples_to_biological_conditions_path analysis_results_directory run_logs_directory
```

For more information and options run the following:
```bash
python3 IgOmeProfiling_pipeline.py -h
```

The code assumes existence of ```bash``` in order to run.  
This can be changed via ```local_command_prefix``` variable in ```global_params.py```.

## Testing
The code contains mock data for testing.

Run the following:
```bash
# From project directory
mkdir output && mkdir output/analysis && mkdir output/logs
python3 IgOmeProfiling_pipeline.py mock_data/exp12_10M_rows.fastq.gz mock_data/barcode2samplename.txt mock_data/samplename2biologicalcondition.txt output/analysis output/logs
```

The entire pipeline might run a few minutes? (TBD).  
The output of the pipeline is (TBD)

## Docker
The code can be containerized using Docker:
```bash
# From project directory
./build_docker.sh
```

The image is built without data and command argument of "-h".
To run the mock-data using docker:
```bash
# From project directory
mkdir test && mkdir test/analysis && mkdir test/logs
docker run --name igome --rm -v ./some_data:/data -v ./test:/output webiks/igome-profile /data/exp12_10M_rows.fastq.gz /data/barcode2samplename.txt /data/samplename2biologicalcondition.txt /output/analysis /output/logs
```

The mock data is included in the docker, to run it:
```bash
mkdir test && mkdir test/analysis && mkdir test/logs
docker run --name igome --rm -v ./test:/output webiks/igome-profile ./mock_data/exp12_10M_rows.fastq.gz ./mock_data/barcode2samplename.txt ./mock_data/samplename2biologicalcondition.txt /output/analysis /output/logs
```

Upload to AWS (using aws-cli with credentials set):
```bash
./aws_upload.sh
```

In AWS machine (with aws-cli credentials set):
```bash
$(aws ecr get-login --no-include-email --region us-west-2)
docker pull 223455578796.dkr.ecr.us-west-2.amazonaws.com/igome-profile:latest
mkdir test && mkdir test/analysis && mkdir test/logs
docker run --name igome --rm -v ./test:/output 223455578796.dkr.ecr.us-west-2.amazonaws.com/igome-profile:latest ./mock_data/exp12_10M_rows.fastq.gz ./mock_data/barcode2samplename.txt ./mock_data/samplename2biologicalcondition.txt /output/analysis /output/logs
```

## Running on multiple machines
TODO allow to set run_using_celery and run_local_in_parallel_mode using env vars
TODO: requirement shared mount at same path, using celery+configuration via env vars, rabbitmq (show how to setup using docker+monitoring using flower). how to run worker (+run using docker)
docker run -d --name rabbit -p 5672:5672 -p 15672:15672 rabbitmq:management
docker run -d --link rabbit:rabbit --name flower -p 5555:5555 mher/flower flower --broker=pyamqp://guest@rabbit//

docker run -d --link rabbit:rabbit -e CELERY_BROKER_URL="pyamqp://guest@rabbit//" -v ./test3:/test3 --name igomeworker --rm webiks/igome-profile-worker
docker run --name igome --rm --link rabbit:rabbit -e CELERY_BROKER_URL="pyamqp://guest@rabbit//" -v ./test3:/test3 webiks/igome-profile:latest ./mock_data/exp12_10M_rows.fastq.gz ./mock_data/barcode2samplename.txt ./mock_data/samplename2biologicalcondition.txt /test3/analysis /test3/logs

## License
TBD