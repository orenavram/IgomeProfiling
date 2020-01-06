#ifndef ___NEEDLEMAN_WUNSCH_H
#define ___NEEDLEMAN_WUNSCH_H

#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>
#include <climits>
#include "MSA.h"
#include "sim_matrices.h"
using namespace std;

class needleman_wunsch
{
public:

	explicit needleman_wunsch(const vector<string> & seqArray,
		int similarity_mode, // when this is 1, we align nucleotide. When 2, AA using Blosum62.
		int match_score = 2,
		int mismatch_score = -1,
		int gap_open = 5,
		int gap_extend = 1);

	explicit needleman_wunsch(//empty constructor, without sequences.
		int similarity_mode, // when this is 1, we align nucleotide. When 2, AA using Blosum62.
		int match_score = 2,
		int mismatch_score = -1,
		int gap_open = 5,
		int gap_extend = 1);

	int getNumberOfSequences() const {return _numberOfSequences;} 
	~needleman_wunsch() {}
	
	void setGapOpenScore(int x) { _gap_open = x; }
	void setGapExtensionScore(int x) { _gap_extend = x; }
	void computeAllPairsMSAs(vector<MSA> & pairwiseMsasToFill);
	MSA computeAffineNWForPair(const string & A, const string & B); // if gap_open = gap_extend --> same as NW without affine gap penalty
	
private:
	vector<string> _seqs;
	sim_matrices _simMat;
	size_t _numberOfSequences;
	int _similarity_mode;
	int _match_score;
	int _mismatch_score;
	int _gap_open;
	int _gap_extend;
};

#endif