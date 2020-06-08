#include "read_PSSMs_from_MAST_file.h"
#include <sstream>
#include <iterator>     // std::istream_iterator
#include <regex>
#include <fstream>

using namespace std;
void fromFileToVectorOfString(const string fileName, vector<string> & allLines);

readPSSM_info_from_file::readPSSM_info_from_file(string pssmFileName) {
	vector<string> allLines;
	fromFileToVectorOfString(pssmFileName, allLines);
	get_aaFreq_and_CorrChar(allLines);
	readFileToPSSM_array(allLines);
	calculate_consensus_seq_for_PSSM_and_max_match_score();
};

void readPSSM_info_from_file::get_aaFreq_and_CorrChar(vector<string> & allLines)
{
	int num_motifs = 0;
	for (size_t i = 0; i<allLines.size(); ++i) {
		string currLine = allLines[i];
		size_t found = currLine.find("Background letter frequencies");
		if (found != string::npos)
		{
			currLine = allLines[i + 1]; // get next line
			stringstream line_element(currLine); // split line by spaces
			istream_iterator<string> begin(line_element);
			istream_iterator<string> end;

			//e.g.:
			//A 0.0625 C 0.03125 D 0.03125 E 0.03125 F 0.03125 G 0.0625 H 0.03125 I 0.03125 K 0.03125 L 0.09375 M 0.03125 N 0.03125 P 0.0625 Q 0.0625 R 0.09375 S 0.09375 T 0.0625 V 0.0625 W 0.03125 Y 0.03125
			vector<string> vstrings(begin, end);
			for (size_t j = 0; j<vstrings.size(); j += 2)
			{
				_alph._correspondingCharacters.push_back(vstrings[j][0]);
				double lol = atof(vstrings[j + 1].c_str());
				_alph._aaFreq.push_back(lol);
			}
			i = allLines.size(); // finish to extract relevant data. Same as break.
		}
	}
	_alph.generateMap();
}

void readPSSM_info_from_file::readFileToPSSM_array(vector<string> & allLines) {
	int num_motifs = 0;
	for (size_t i = 0; i<allLines.size(); ++i) {
		//if MOTIF - start inner loop
		string currLine = allLines[i];
		std::size_t found = currLine.find("MOTIF");
		if (found != std::string::npos)
		{
			//std::cout << "found MOTIF at line " << i << "this is motif number " << num_motifs << endl;
			string motif_name = currLine.substr(6); // exclude from the name the beginning "MOTIF ".
			PSSM crrPSSM(_alph);
			crrPSSM.setName(motif_name);
			// read _Nsite
			currLine = allLines[i + 1];
			smatch match;
			regex rgx(".*nsites=\\s+([\\d]+)");
			if (regex_search(currLine, match, rgx))
			{
				string Nsite_str = match[1].str();
				crrPSSM.set_Nsite(stod(Nsite_str)); // Nsite is the number of unique peptides inducing the motif.
			}
			i = i + 2; // let's collect matrix lines
			vector<string> PSSM_Block;
			currLine = allLines[i];
			while (currLine.compare("") != 0) // did not reach empty line - end of block
			{
				//std::cout << currLine << endl;
				PSSM_Block.push_back(currLine);
				++i;
				if (i<allLines.size())
				{
					currLine = allLines[i];
				}
			}
			crrPSSM.setMatrix(PSSM_Block);
			crrPSSM.Add_PseudoCount(1);
			_PSSM_array.push_back(crrPSSM);
			pair<map<string, int>::iterator, bool> ret;
			ret = _PSSM_Name_To_ArrayIndex.insert(std::pair<string, int>(motif_name, num_motifs));
			if (ret.second == false) {
				cout << "Motif name: '" << motif_name << "' already existed ";
				cout << " with a value of " << ret.first->second << "please make sure the motifs names are unique!!.." << endl;;
			}
			//std::cout << "Motif is '" << crrPSSM.PSSM_name << "'" << endl;
			++num_motifs;
		}
	}
}

void readPSSM_info_from_file::calculate_consensus_seq_for_PSSM_and_max_match_score() {
	for (size_t i = 0; i < _PSSM_array.size(); ++i)
	{
		vector<size_t> consensus = _PSSM_array[i].GetConsensusSeq();
		int best_match_start = 0;
		string best_match_string = "";
		double consensus_seq_score = _PSSM_array[i].computeScore(consensus, best_match_start);
		_PSSM_array[i].setConsensus(consensus);
		_PSSM_array[i].setMaxScore(consensus_seq_score);
	}
}

void readPSSM_info_from_file::update_PSSM_cutoff_from_file (const string CutoffsFileName) {
	if (!ifstream(CutoffsFileName))
	{
		std::cout << "The cutoffs Per PSSM file named: " << CutoffsFileName << "does not exists!\n you need to run the program with the -calc_cutoff flag to calculate it..." <<endl;
		exit(4);
	}
	ifstream file(CutoffsFileName);
	string   line;
	getline(file, line); // read line by line
	while (!file.eof()) // file not finished
	{
		if (line.substr(0, 4).compare("###\t") == 0) // got to new motif name
		{
			stringstream   linestream(line); // create a stream
			string tmp;
			string MotifName;
			getline(linestream, tmp, '\t');
			getline(linestream, MotifName, '\t');
			while ((getline(file, line)) && (line.substr(0, 4).compare("###\t") != 0)) // read all data lines
			{
				stringstream   linestream(line); // create a stream
				string         data;
				string SeqType;
				string cutoff;
				getline(linestream, SeqType, '\t'); // read up-to the first tab (discard tab).
				getline(linestream, cutoff, '\t');
				// update the relevant PSSM cutoff
				map<string, int>::iterator p;
				p = _PSSM_Name_To_ArrayIndex.find(MotifName);
				int PSSM_index;
				if (p != _PSSM_Name_To_ArrayIndex.end())
				{
					PSSM_index = p->second;
					_PSSM_array[PSSM_index].set_SeqTypCutoff(SeqType, stod(cutoff));
				}
				else
				{
					cout << "Can't find the motif named '" << MotifName << "' in the PSSM array... ignotr it cutoff score from the input file!!" << endl;
				}

			}
		}
	}
}
