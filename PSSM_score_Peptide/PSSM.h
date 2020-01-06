#ifndef PSSM_____H
#define PSSM_____H

#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <map>
#include <algorithm>    // std::min_element, std::max_element
#include <iostream>
#include <cmath>  // because of unix

using namespace std;
#include "alphabet.h"

class PSSM {
public:
	
	
	explicit PSSM(alphabet& alph) :_alph(alph){set_Nsite(0);}
	~PSSM(){};
	PSSM(const PSSM & p1) :_alph(p1._alph) { PSSMmatrix = p1.PSSMmatrix;
							PSSM_name  = p1.PSSM_name;
							_PSSM_consensus_seq = p1._PSSM_consensus_seq;
							PSSM_MaxScore = p1.PSSM_MaxScore;
							_Nsite = p1._Nsite;
							dist_of_random_peptides_scores.clear();
							seq_type_cutoff = p1.seq_type_cutoff;
							}

	void setName( string name ){ PSSM_name = name;}
	void setMatrix(const vector<string> & PSSMLines);
	void setConsensus (const  vector<size_t> &consensus_seq){ _PSSM_consensus_seq =consensus_seq;}
	void setMaxScore (const double score){PSSM_MaxScore=score;}
	void set_SeqTypCutoff(const string SeqType,const double cutoff){seq_type_cutoff[SeqType]=cutoff;}
	vector<size_t> GetConsensusSeq() const;
	void set_Nsite (double Nsite){_Nsite=Nsite;}
	void Add_PseudoCount(double PseudoCountSize);
	double computeScore(const vector <size_t> & seq, int & best_match_start) const;
	double computeScoreExactr(const vector<size_t>& seq_string, size_t startPos, size_t endPos) const;
	double computeScoreExacrPos(size_t i, const size_t c) const;
	size_t PSSM_length() const {
		return PSSMmatrix.size();
	
	}
	PSSM randomize();
	//int integerOfChar(const char s) const;

	//members:
	vector< vector<double> > PSSMmatrix;
	string PSSM_name;
	vector<size_t> _PSSM_consensus_seq;
	double PSSM_MaxScore;
	double _Nsite;
	map<string, vector<double>> dist_of_random_peptides_scores;
	map<string, double> seq_type_cutoff;
	alphabet& _alph;

};

#endif
