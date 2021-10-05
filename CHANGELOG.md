# Change Log
All notable changes to this project will be documented in this file.

### [1.2.0] - 2021-10-05

### Changes
- Changed defualt values of parameters:
    - Phase 1: 
        - Changed maximum_length_required from 12 to 14
        - Calculate rpm by defualt unless using flag no_calculate_rpm
    - Phase 2:
        - Changed aln_cutoff from 20 to 24
        - Changed pcc_cutoff from 0.6 to 0.7
        - Changed threshold from 0.5 to 0.6
        - Changed word_length from 2 to 4
        - Changed discard from 1 to 4
    - Phase 3:    
        - Changed shuffles from 5 to 10
        - Scanning with using rpm faa file unless using no_use_rpm_faa_scanning
        - Log the hit sequences while scanning unless using no_output_sequences_scanning
- Changed the position of stop machines AWS to be after create done file. 

### Added
- Added CHANGELOG file.

### [1.1.0] - 2021-09-22

### Added
- Added new file - positive motifs: keep only positive motifs before the random forest.

### [1.0.0] - 2021-08-01

- 🎉 first stable release!
 
