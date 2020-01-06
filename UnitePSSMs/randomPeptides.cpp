#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <random> // random number 
#include <fstream> 
#include <string>
#include <algorithm>
using namespace std;

#include "randomPeptides.h"
randomPeptide::randomPeptide(vector<double> characterFrequencies, vector<char> correspondingCharacters, size_t numberOfSeqToSimulate, size_t sequenceLength) : _characterFrequencies(characterFrequencies), _correspondingCharacters(correspondingCharacters), _numberOfSeqToSimulate(numberOfSeqToSimulate), _sequenceLength(sequenceLength),_CysLoop (0) { // deafualt constroctor - NoCysLoop
     srand(static_cast<unsigned int>(time(NULL)));
}
	
randomPeptide::randomPeptide(vector<double> characterFrequencies, vector<char> correspondingCharacters, size_t numberOfSeqToSimulate, size_t sequenceLength,bool CysLoop) : _characterFrequencies(characterFrequencies), _correspondingCharacters(correspondingCharacters), _numberOfSeqToSimulate(numberOfSeqToSimulate), _sequenceLength(sequenceLength),_CysLoop (CysLoop) {
     srand(static_cast<unsigned int>(time(NULL)));
}	
void randomPeptide::generateRandomSequences() {
//	std::mt19937 Seq_eng;  // a core engine class
//	std::random_device dev_random;
//	Seq_eng.seed(dev_random()); //real random seed
//	std::uniform_real_distribution<double> unif; //std::uniform_real<double> unif(0, 1);
	for (size_t i=0; i<_numberOfSeqToSimulate; ++i) {
		string tmp;
		if (_CysLoop)
		{
			tmp="C";
		}
		tmp.append(simulateOneSeq());
		if (_CysLoop)
		{
			tmp.append("C");
		}
		_simulatedSequences.push_back(tmp);
	}
}

string randomPeptide::simulateOneSeq() {
	string tmp;
	for (size_t i=0; i<_sequenceLength; ++i) {
		char c = simulateOneCharacter();
		tmp += c;
	}
	return tmp;
}

char randomPeptide::simulateOneCharacter() {
	int Rand_int = rand();
	double Rand =((double) Rand_int / (RAND_MAX));
	int i=0;
	while((i<_characterFrequencies.size()) && ( Rand >= _characterFrequencies[i]))
	{
		Rand=Rand-_characterFrequencies[i];
		++i;
//		if (i>19)
//		{
//			cout<<"OUT OF RANGE"<<endl;
//		}
	}
	if (i==_characterFrequencies.size()){--i;} // to avoid out of range when the Rand is getting to 0
	return _correspondingCharacters[i];

}

