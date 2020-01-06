#ifndef PSSM_____H
#define PSSM_____H

#include <vector>
#include <string>


using namespace std;


class PSSM {
public:	
	PSSM(size_t numberOfSupportingSequences, const string& pssmName, const vector<string> & PSSMLines, const vector<char> & correspondingCharacters, const vector<double> & aaFreq, size_t pseudoCountValue = 1); //getting the PSSM from a section of an input file.
	~PSSM(){};
	void setName(string name ){ _PSSM_name = name;}
	string computeConsensusSeq(); // computes it,  returns it, and keeps it in _PSSM_consensus_seq
	string getConsensus() const {
		return _PSSM_consensus_seq;
	}
	void set_Nsite (size_t Nsite){_Nsite=Nsite;}
	size_t get_Nsite() const { return _Nsite; };
	string get_PSSM_name() const { return _PSSM_name; };
	void addPseudoCount(size_t val); // adding val to each entry in the PSSM.
	
	double getValue(size_t position, size_t column) const {
		return _PSSMmatrix[position][column];
	}

	double getLogValue (size_t position, size_t column) const {
		return _PSSMmatrixLOG[position][column];
	}

	size_t alphabetSize() const {
		if (_PSSMmatrix.size() == 0) return 0;
		return _PSSMmatrix[0].size();
	}

	// This function returns the percentile of the score. 
	// If percentile is 0.95, it will give the score X for which 95% of random peptides have score less than X.
	// As often the peptide size is less than the PSSM, this score also depends on the peptide size.
	// For example, when the PSSM is of length 40, and the peptide size is of length 5, we try all windowns of
	// size 5 for each peptide.
	// As this procedure relies on simulations, we also provide the number of simulations needed to estimate this percentile.
	double giveScorePercentile(size_t peptideSize, double percentile, size_t numberOfSimulations);

private:
	vector< vector<double> > _PSSMmatrix; // positions times alphabet
	vector< vector<double> > _PSSMmatrixLOG; // positions times alphabet
	string _PSSM_name;
	string _PSSM_consensus_seq;
	double _PSSM_MaxScore; // what is the max possible score of the PSSM (the score of the consensus seq).
	size_t _Nsite; // how many sequences were used to build the PSSM
	vector<char> _correspondingCharacters; // for a given column in the PSSM, the first row is the probability of the first character, etc.
	vector<double> _aaFreq; // the aa freq in the background model.
	vector<double> _cummulative_aaFreq; // the aa freq in the background model.
	double getScoreOf_a_RandomPeptide(const vector<size_t> & randomPeptideEntries, size_t motifSize); // internal procedure to compute the score of a random peptide

				   //map<string,vector<double>> dist_of_random_peptides_scores; 
				   //map<string,double> seq_type_cutoff;
	

};

#endif
