#ifndef RANDOM_PEPTIDES_____H
#define RANDOM_PEPTIDES_____H

class randomPeptides {
public:
	//randomPeptides(vector<double> characterFrequencies, vector<char> correspondingCharacters, size_t numberOfSeqToSimulate, size_t sequenceLength);
	//randomPeptides(vector<double> characterFrequencies, vector<char> correspondingCharacters, size_t numberOfSeqToSimulate, size_t sequenceLength,bool CysLoop);
	randomPeptides(vector<double>& characterFrequencies, size_t numberOfSeqToSimulate, size_t sequenceLength);
	randomPeptides(vector<double>& characterFrequencies, size_t numberOfSeqToSimulate, size_t sequenceLength,bool CysLoop);
	void generateRandomSequences();
	void simulateOneSeq(vector<size_t> &tmp);
	size_t simulateOneCharacter();

	vector<double> _characterFrequencies;
	//vector<char> _correspondingCharacters;
	size_t _numberOfSeqToSimulate;
	size_t _sequenceLength;
	vector<vector<size_t> > _simulatedSequences;
	bool _CysLoop;
};
#endif