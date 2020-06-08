#ifndef READ_PSSM_INFO_FROM_FILE_H
#define READ_PSSM_INFO_FROM_FILE_H

#include <string>
#include <vector>
#include <map>
using namespace std;

#include "PSSM.h"
#include "alphabet.h"

class readPSSM_info_from_file {
public:
	readPSSM_info_from_file(string pssmFileName);
	void get_aaFreq_and_CorrChar(vector<string> & allLines);
	void readFileToPSSM_array(vector<string> & allLines);
	void calculate_consensus_seq_for_PSSM_and_max_match_score();
	void update_PSSM_cutoff_from_file(string CutoffsFileName);



	alphabet _alph;
	vector<PSSM> _PSSM_array;
	map<string, int> _PSSM_Name_To_ArrayIndex; // A map between the motif name and its index.
};



#endif
