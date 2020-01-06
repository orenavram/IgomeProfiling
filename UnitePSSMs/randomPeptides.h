#ifndef RANDOM_PEPTIDES_____H
#define RANDOM_PEPTIDES_____H

class randomPeptide {
public:
	randomPeptide(vector<double> characterFrequencies, vector<char> correspondingCharacters, size_t numberOfSeqToSimulate, size_t sequenceLength);
	randomPeptide(vector<double> characterFrequencies, vector<char> correspondingCharacters, size_t numberOfSeqToSimulate, size_t sequenceLength,bool CysLoop);
	void generateRandomSequences();
	string simulateOneSeq();
	char simulateOneCharacter();

	vector<double> _characterFrequencies;
	vector<char> _correspondingCharacters;
	size_t _numberOfSeqToSimulate;
	size_t _sequenceLength;
	vector<string> _simulatedSequences;
	bool _CysLoop;
};
#endif