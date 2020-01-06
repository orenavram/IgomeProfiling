#include "fromFileToPSSM_array.h"

#include <sstream>
#include <fstream>
#include <regex>
#include <iostream>
using namespace std;

void fromFileToVectorOfString(const string fileName, vector<string> & allLines) {
	ifstream myfile;
	myfile.open(fileName.c_str());
	if (myfile.fail()) { 
		cerr << "Failed to open PSSMs file " << fileName << endl;
		cout << "Failed to open PSSMs file " << fileName << endl;
		exit(EXIT_FAILURE);
	}
	string tmp;
	while (getline(myfile, tmp)) {
		allLines.push_back(tmp);
	}
	myfile.close();
}

void readFileToPSSM_array(const string fileName, vector<PSSM> & PSSM_array, map<string, int> & PSSM_Name_To_ArrayIndex) {
	// TO DO: CREATE ALSO A MAP BETWEEN PSSM NAME AND ARRAY INDEX IN PSSM_ARRAY...
	vector<string> allLines;
	fromFileToVectorOfString(fileName, allLines);


	// read aa frequenices and character array:
	vector <double> aaFreq;
	vector <char> correspondingCharacters;
	int num_motifs = 0;
	for (size_t i = 0; i<allLines.size(); ++i) {
		string currLine = allLines[i];
		std::size_t found = currLine.find("Background letter frequencies");
		if (found != std::string::npos)
		{
			currLine = allLines[i + 1]; // get next line
			stringstream line_element(currLine); // split line by spaces
			istream_iterator<string> begin(line_element);
			istream_iterator<string> end;
			vector<string> vstrings(begin, end);
			for (size_t j = 0; j<vstrings.size(); j += 2)
			{
				correspondingCharacters.push_back(vstrings[j][0]);
				double lol = atof(vstrings[j + 1].c_str());
				aaFreq.push_back(lol);
			}
			i = allLines.size(); // finish to extract relevant data
		}
	}

	// here we read ther main PSSM data
	num_motifs = 0;
	for (size_t i = 0; i<allLines.size(); ++i) {
		//if MOTIF - start inner loop
		string currLine = allLines[i];
		std::size_t found = currLine.find("MOTIF");
		if (found != std::string::npos)
		{
			//std::cout << "found MOTIF at line " << i << "this is motif number " << num_motifs << endl;
			size_t numberOfSites;
			string motif_name = currLine.substr(6);
			
			// read _Nsite
			currLine = allLines[i + 1];
			smatch match;
			regex rgx(".*nsites=\\s+([\\d]+)");
			if (regex_search(currLine, match, rgx))
			{
				string Nsite_str = match[1].str();
				numberOfSites = stoi(Nsite_str);
			}
			i = i + 2; // let's collect matrix lines
			vector <string> PSSM_Block;
			currLine = allLines[i];
//			cout << "**************************************" << endl;
			while (currLine.compare("") != 0) // did not reach empty line - end of block
			{
//				std::cout << currLine << endl;
				PSSM_Block.push_back(currLine);
				++i;
				if (i<allLines.size())
				{
					currLine = allLines[i];
				}
			}
//			std::cout << "**************************************" << endl;
			PSSM crrPSSM(numberOfSites, motif_name, PSSM_Block,correspondingCharacters, aaFreq);
			
			PSSM_array.push_back(crrPSSM);
			std::pair<std::map<string, int>::iterator, bool> ret;
			ret = PSSM_Name_To_ArrayIndex.insert(std::pair<string, int>(motif_name, num_motifs));
			if (ret.second == false) {
//				std::cout << "Motif name: '" << motif_name << "' already existed ";
//				std::cout << " with a value of " << ret.first->second << "please make sure the motifs names are unique!!.." << endl;;
			}
			//std::cout << "Motif is '" << crrPSSM.PSSM_name << "'" << endl;
			++num_motifs;
		}
	}
}