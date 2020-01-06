#ifndef ALPHABET_____H
#define ALPHABET_____H

#include <vector>
#include <map>
using namespace std;

#define GAP 20

class alphabet {
public:
	alphabet(vector<char> & correspondingCharacters, vector<double>& aaFreq) : _correspondingCharacters(correspondingCharacters), _aaFreq(aaFreq) {
		generateMap();
	}
	alphabet(){} // empty constructor

	void generateMap() {
		for (size_t k = 0; k < _correspondingCharacters.size(); ++k) {
			_alphabetMap[_correspondingCharacters[k]] = k;
		}
		_alphabetMap['-'] = GAP;
	}
	vector<char> _correspondingCharacters;
	vector<double> _aaFreq;
	map<char, size_t> _alphabetMap;
};

#endif
