# Cross experiments IgOme Profiling pipeline

### Phase 1 - reads filtration
Run multi experiments of phase 1 by adding json file with the flag "--multi_exp_config_reads".

Structure of the json file:
```
  {
      "configuration":{},
      "runs":{
          "name_exp":{
              "reads_path": "",
              "fastq": "",
              "barcode2sample": "",
              "done_path": ""
          }
      }
  }
```
The four parameters under name_exp are mast.

### Phase 2 - motif inference
Run multi experiments of phase 2 by adding json file with the flag "--multi_exp_config_inference".

Structure of the json file:
```
  {
      "configuration":{},
      "runs":{
          "name_exp":{
              "reads_path": "",
              "motifs_path": "",
              "sample2bc": "",
              "done_path": ""
          }
      }
  }
```
The four parameters under name_exp are mast.

### Phase 3 - model fitting
Run multi experiments and create cross between multi experiments by adding json file with the flag "--cross_experiments_config".

Structure of the json file:
```
  {
      "configuration":{},
      "runs":{
          "name_exp":{
              "configuration": {
                "reads_path": "",
                "motifs_path": "",
                "sample2bc": "",
                "done_path": ""
              },
              "cross": {
                  "name_cross": {
                        "motifs": [],
                        "samples": {
                            "BC_name": {
                                "biologicalConditions":[]
                            },
                            "other": {
                                "biologicalConditions":[]
                            }
                        }
                  }
              }
          }
      }
  }
```
The four parameters under name_exp/confiduration are mast.
