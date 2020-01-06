#include "MSA.h"

void MSA::unalignSeqs()
{
	for(size_t j=0; j<_numberOfSequences; j++)	
	{
		string curr_unaligned = "";
		for(size_t i=0; i<getMSAlength(); i++) 
		{
			if(_alignedSeqs[j][i]!='-')
			{
				curr_unaligned += _alignedSeqs[j][i];
			}
		}
		_unalignedSeqs.push_back(curr_unaligned);
	}
}

MSA::MSA(const MSA& other )//copy constructor
{
	_alignedSeqs = other._alignedSeqs;
	_unalignedSeqs = other._unalignedSeqs;
	_numberOfSequences = other._numberOfSequences;
	_score = other._score;
}

MSA&  MSA::operator=(const MSA& other)
{
    _alignedSeqs= other._alignedSeqs;
	_unalignedSeqs = other._unalignedSeqs;
	_numberOfSequences=other._numberOfSequences;
	_score = other._score;
	return *this;
}


bool MSA::isMatch(size_t pos) const {
	bool match = true;
	for (size_t i = 0; i < _alignedSeqs.size(); ++i) {
		if (_alignedSeqs[i].at(pos) == '-') {
			match = false;
			break;
		}
	}
	return match;
}

MSA::~MSA()
{
	_alignedSeqs.clear();
}
