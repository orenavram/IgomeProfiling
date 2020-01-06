#include "PSSM.h"
#include "SEQ.h"

#include <algorithm>
using namespace std;

void PSSM::setMatrix(const vector<string> & PSSMLines)
{
	size_t numberOfPositions=PSSMLines.size(); // init lines of matrix
	PSSMmatrix.resize(numberOfPositions);
	for (size_t MotifPos=0;MotifPos<PSSMLines.size();MotifPos++) // fill in lines
	{
		string line=PSSMLines[MotifPos]; //s = "0 0.926111458985598 0.000626174076393237 0 0.00438321853475266 0.0050093926111459 0 0.00125234815278647 0.00187852222917971 0.00250469630557295 0.00250469630557295 0 0.000626174076393237 0 0.00438321853475266 0.0219160926737633 0.00751408891671885 0 0.0137758296806512 0.00751408891671885";
		stringstream line_element(line); // split line by spaces
		istream_iterator<string> begin(line_element);
		istream_iterator<string> end;
		vector<string> vstrings(begin, end);
		PSSMmatrix[MotifPos].resize(vstrings.size()+1); // last position is for padding gap character
		for (size_t i=0;i<vstrings.size();i++)
		{
			double lol = atof(vstrings[i].c_str()); 
			PSSMmatrix[MotifPos][i]=lol;
		}
		PSSMmatrix[MotifPos][vstrings.size()]=0; // last position is for padding gap character
	}

	//cout<<"done with pssm matrix before pseudocounts"<<endl;
	//cout<<"num elements in first row for example is "<<PSSMmatrix[0].size()<<endl;
	//cout<<"first row is:"<<endl;
	double sum_first_row = 0;
	for(size_t i = 0; i < PSSMmatrix[0].size(); i++)
	{
		sum_first_row = sum_first_row + PSSMmatrix[0][i];
		//cout<<PSSMmatrix[0][i]<<endl;
	}
	//cout<<"sum of elements in first row for example is "<<sum_first_row<<endl;
	//cout<<"gap element in first row for example is "<<PSSMmatrix[0][PSSMmatrix[0].size()-1]<<endl;


}

vector<size_t> PSSM::GetConsensusSeq() const {
	vector<size_t> ConsensusSeq;
	for (size_t MotifPos=0;MotifPos<PSSMmatrix.size();MotifPos++) {
		size_t MaxCharIndex=distance(PSSMmatrix[MotifPos].begin(), max_element(PSSMmatrix[MotifPos].begin(),PSSMmatrix[MotifPos].end()));
		ConsensusSeq.push_back(MaxCharIndex);
	}
	return ConsensusSeq;
}

void PSSM::Add_PseudoCount(double PseudoCountSize)
{
	if (_Nsite==0)
	{
		cout<<"ERROR: Cna't employ AddPseudoCount on motif: '"<<PSSM_name<<"' because Nsite is not defiend..."<<endl;
		return;
	}
	else
	{
		for (size_t MotifPos=0;MotifPos<PSSMmatrix.size();MotifPos++) 
		{
			for (size_t i=0;i<PSSMmatrix[MotifPos].size();i++)
			{
				double currValue=PSSMmatrix[MotifPos][i];
				double newVal=((currValue*_Nsite)+PseudoCountSize)/(_Nsite+(PseudoCountSize*PSSMmatrix[MotifPos].size()));
//				newVal=log(newVal)/log(2.0);
				PSSMmatrix[MotifPos][i]=newVal;
			}
		}
	}


	//cout<<"num seqs was: "<<_Nsite<<endl;
	//cout<<"done with pssm matrix after pseudocounts"<<endl;
	//cout<<"num elements in first row for example is "<<PSSMmatrix[0].size()<<endl;
	//cout<<"first row is:"<<endl;
	double sum_first_row = 0;
	for(size_t i = 0; i < PSSMmatrix[0].size(); i++)
	{
		sum_first_row = sum_first_row + PSSMmatrix[0][i];
		//cout<<PSSMmatrix[0][i]<<endl;
	}
	//cout<<"sum of elements in first row for example is "<<sum_first_row<<endl;
	//cout<<"gap element in first row for example is "<<PSSMmatrix[0][PSSMmatrix[0].size()-1]<<endl;
}

vector<size_t> PadSeq(const vector<size_t>& seq, size_t PaddingLength); // implemented in the main file.

double PSSM::computeScore(const vector<size_t> &seq, int & best_match_start) const {
	double maxScore = 1; // this exceeds the maximum possible score (sum of log prob values can be max 0)

	best_match_start = 0; // init
	vector<size_t> best_match_string;
	int PSSM_Size = PSSMmatrix.size();
	int paddingSize = PSSM_Size;
	vector<size_t > SeqPadded = PadSeq(seq, paddingSize);
	for (size_t posInSeq = 0; posInSeq < SeqPadded.size() - PSSM_Size + 1; ++posInSeq) {
		//vector<size_t> stringToCheck = SeqPadded.substr(i, PSSM_Size);
		double score = computeScoreExactr(SeqPadded, posInSeq, posInSeq +PSSM_Size);
		if (maxScore == 1) //for first iteration
		{
			maxScore = score;
			//best_match_string = stringToCheck;
		}
		if (score>maxScore)
		{
			maxScore = score;
			// best_match_start = posInSeq;
			best_match_start = posInSeq - paddingSize; // without the padding
			//best_match_string = stringToCheck;
		}
	}
	if (maxScore == 0)
	{
		SEQ tmp(_PSSM_consensus_seq, "", 1, _alph);
		//cout << PSSM_name << "\t" << tmp.getStringOfSeq() << "\t" << PSSMmatrix.size() << "\t" << seq.size() << "\t";
		SEQ tmpSeq(seq,"",0,_alph);
		//cout << tmpSeq.getStringOfSeq();
		//cout<<"\t" << maxScore << endl; // QA
	}
	return maxScore;
}

double PSSM::computeScoreExactr(const vector<size_t>& seq_string, size_t startPosInSeq,size_t endPosInSeq) const {
	double sum = 0;
	for (size_t posInPSSM = 0; posInPSSM < PSSM_length() ; ++posInPSSM) {
		double pos_score = computeScoreExacrPos(posInPSSM,seq_string[startPosInSeq+ posInPSSM]);
		sum += pos_score;
	}
	if (sum == 0) {
		cout << "we have a problem!";
	}
	return sum;
}

double PSSM::computeScoreExacrPos(size_t posInPSSM, const size_t charInSeq) const {
	double res = PSSMmatrix[posInPSSM][charInSeq];
	return log(res);
}

/*int PSSM::integerOfChar(const char s) const {
	switch (s) {
	case 'A': case'a': return 0; break;
	case 'C': case'c': return 1; break;
	case 'D': case'd': return 2; break;
	case 'E': case'e': return 3; break;
	case 'F': case'f': return 4; break;
	case 'G': case'g': return 5; break;
	case 'H': case'h': return 6; break;
	case 'I': case'i': return 7; break;
	case 'K': case'k': return 8; break;
	case 'L': case'l': return 9; break;
	case 'M': case'm': return 10; break;
	case 'N': case'n': return 11; break;
	case 'P': case'p': return 12; break;
	case 'Q': case'q': return 13; break;
	case 'R': case'r': return 14; break;
	case 'S': case's': return 15; break;
	case 'T': case't': return 16; break;
	case 'V': case'v': return 17; break;
	case 'W': case'w': return 18; break;
	case 'Y': case'y': return 19; break;

	case '-': case'_': return 20; break; // padding character

	default:
		vector<string> err;
		err.push_back(" The amino-acid sequences contained the character: ");
		err[0] += s;
		err.push_back(" Amino acid was not one of the following: ");
		err.push_back(" A, B, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, -");
		err.push_back(" a, b, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v, _");
		for (size_t k = 0; k < err.size(); ++k) {
			cerr << err[k];
		}
		cerr << endl;
	}// end of switch
	return -99; // never supposed to be here.	
}// end of function
*/

PSSM PSSM::randomize() {
	PSSM res(*this);
	random_shuffle(res.PSSMmatrix.begin(), res.PSSMmatrix.end());
	return res;
}