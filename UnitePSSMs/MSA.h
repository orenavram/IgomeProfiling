#ifndef _MSA
#define _MSA
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <map>
#include <iostream> // for printing the MSA to cout.
using namespace std;

class MSA
{
public:

	MSA(const vector<string> & seqArray,const int score): _alignedSeqs(seqArray), _score(score)
	{
		_numberOfSequences = _alignedSeqs.size();
		unalignSeqs();
	} //constructor from vector string
	MSA(){};
	MSA(const MSA& other);   //copy constructor

	size_t getMSAlength() const {return _alignedSeqs[0].size();}
	size_t getNumberOfSequences() const {return _numberOfSequences;}
	int getAlignmentScore() const { return _score; }

	vector<string> getUnalignedSeqs () const {return _unalignedSeqs;}
	vector<string> getAlignedSeqs () const {return _alignedSeqs;}
	void printAligned() {
		for (size_t i = 0; i < _alignedSeqs.size(); ++i) {
			cout << _alignedSeqs[i] << endl;
		}
	}
	
	~MSA();
	MSA& operator=(const MSA& );

	bool isMatch(size_t pos) const;
	bool isGap(size_t sequenceNumber, size_t pos) const {
		return (_alignedSeqs[sequenceNumber].at(pos) == '-');
	}

private:
	size_t _numberOfSequences; // NUMBER OF SEQUENCES IN THE MSA
	vector<string> _alignedSeqs; //The aligned sequences
	vector<string> _unalignedSeqs; //The unaligned sequences
	int _score;


	//methods
	void unalignSeqs();
};
#endif