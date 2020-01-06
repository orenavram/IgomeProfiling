#include "HIT.h"

void HIT::set_hit_info (double const & hit_score, int const & match_pos)
{
	_match_score = hit_score;
	_match_pos = match_pos;
}
