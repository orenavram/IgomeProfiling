#ifndef HIT_____H
#define HIT_____H

#include <vector>
#include <string>
#include "SEQ.h"
#include <iostream>
#include <sstream>      // std::stringstream
#include <regex>

using namespace std;

class HIT {
public:
	SEQ _seq;
	double _match_score;
	int _match_pos;
	
	HIT(const SEQ & HitSEQ) : _seq(HitSEQ)
	{
		
		_match_score=-1.0;
		_match_pos=-1;
	} //constructor
	~HIT(){};
	void set_hit_info (double const & hit_score, int const & match_pos);
};

#endif