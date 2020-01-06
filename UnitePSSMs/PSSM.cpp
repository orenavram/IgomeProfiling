#include "PSSM.h"
#include <sstream>
#include <iterator>
#include <map>
#include <algorithm>    // std::min_element, std::max_element
#include <iostream>
#include <random>
using namespace std;

string PSSM::computeConsensusSeq()
{
string ConsensusSeq;
for (size_t MotifPos = 0; MotifPos < _PSSMmatrix.size(); MotifPos++) {
	int MaxCharIndex = distance(_PSSMmatrix[MotifPos].begin(), max_element(_PSSMmatrix[MotifPos].begin(), _PSSMmatrix[MotifPos].end()));
	ConsensusSeq += _correspondingCharacters[MaxCharIndex];
}
return ConsensusSeq;
}

// -------------------------------------------------------------
// Here is an example of a possible input
//0 0.914503816793893 0 0.00152671755725191 0.0183206106870229 0.00763358778625954 0.00152671755725191 0 0 0 0.00458015267175573 0 0.00152671755725191 0 0.0137404580152672 0.0122137404580153 0 0 0.0106870229007634 0.0137404580152672
//0.00563380281690141 0 0.0028169014084507 0.190140845070423 0 0.0211267605633803 0.0112676056338028 0.00422535211267606 0.00845070422535211 0.0154929577464789 0.311267605633803 0.00140845070422535 0.359154929577465 0.0154929577464789 0.0112676056338028 0.00845070422535211 0.0112676056338028 0.0169014084507042 0.00140845070422535 0.00422535211267606
//0.0053475935828877 0.0267379679144385 0.00320855614973262 0.00748663101604278 0.00748663101604278 0.00427807486631016 0.00855614973262032 0.0203208556149733 0.00427807486631016 0.281283422459893 0.303743315508021 0.00106951871657754 0.0203208556149733 0.013903743315508 0.00427807486631016 0.027807486631016 0.00962566844919786 0.247058823529412 0 0.00320855614973262
//0 0.0337078651685393 0.0112359550561798 0.00306435137895812 0.00510725229826353 0.00102145045965271 0.00817160367722165 0.0122574055158325 0.0183861082737487 0.438202247191011 0.274770173646578 0.00715015321756895 0.024514811031665 0.00715015321756895 0.0153217568947906 0.0112359550561798 0.0868232890704801 0.0326864147088866 0.00919305413687436 0
//0.00305188199389624 0.00101729399796541 0.00712105798575788 0.00101729399796541 0.0081383519837233 0 0.247202441505595 0.452695829094608 0 0.140386571719227 0.0183112919633774 0.00712105798575788 0.0305188199389624 0.0081383519837233 0.00305188199389624 0.00915564598168871 0.0101729399796541 0.0406917599186165 0.00203458799593082 0.0101729399796541
//0 0.00101626016260163 0 0.0040650406504065 0.0152439024390244 0 0.0345528455284553 0.226626016260163 0.0101626016260163 0.474593495934959 0.016260162601626 0.0111788617886179 0.00711382113821138 0.144308943089431 0.00304878048780488 0.0182926829268293 0.00914634146341463 0.0111788617886179 0.00914634146341463 0.0040650406504065
//0.0101317122593718 0 0.803444782168186 0.0840932117527862 0 0.022289766970618 0.0172239108409321 0.00202634245187437 0.00506585612968592 0.00202634245187437 0.00101317122593718 0.0131712259371834 0 0.0060790273556231 0 0.0162107396149949 0 0.0101317122593718 0 0.00709219858156028
//0 0.00100502512562814 0 0 0.00904522613065327 0 0.00703517587939699 0.0120603015075377 0 0.923618090452261 0.00100502512562814 0 0.0170854271356784 0 0.00804020100502513 0.00402010050251256 0.00201005025125628 0.0130653266331658 0.00201005025125628 0
//0.00802407221664995 0 0.925777331995988 0.0120361083249749 0 0.0170511534603811 0.00702106318956871 0 0.00100300902708124 0 0 0.0100300902708124 0.00401203610832497 0.00100300902708124 0.00100300902708124 0.00200601805416249 0 0.00702106318956871 0 0.00401203610832497
//0 0 0.0101214574898785 0.00101214574898785 0.01417004048583 0 0.922064777327935 0 0 0.00506072874493927 0 0.00809716599190283 0.00708502024291498 0.00708502024291498 0.00809716599190283 0.00303643724696356 0 0.00101214574898785 0 0.0131578947368421
//0.0113168724279835 0.118312757201646 0.00823045267489712 0.00205761316872428 0.00925925925925926 0.0102880658436214 0.00205761316872428 0.0216049382716049 0.00102880658436214 0.0216049382716049 0 0.0195473251028807 0.0154320987654321 0.00205761316872428 0.0123456790123457 0.435185185185185 0.00720164609053498 0.279835390946502 0.0154320987654321 0.00720164609053498
//0.0036231884057971 0.228260869565217 0.0132850241545894 0.0072463768115942 0 0.00845410628019324 0.0108695652173913 0.00120772946859903 0.00483091787439614 0.0144927536231884 0 0.00120772946859903 0.300724637681159 0.309178743961353 0.0205314009661836 0.0205314009661836 0.0144927536231884 0.0108695652173913 0.0181159420289855 0.0120772946859903
//0.0114285714285714 0.0228571428571429 0.00857142857142857 0.0371428571428571 0.00571428571428571 0.0142857142857143 0.0171428571428571 0.00285714285714286 0 0.0314285714285714 0.0114285714285714 0 0.0285714285714286 0.00285714285714286 0.00571428571428571 0.0171428571428571 0.00571428571428571 0.762857142857143 0 0.0142857142857143
//
// Each line has 20 entires and corresponds to one column of the PSSM
PSSM::PSSM(size_t numberOfSupportingSequences, const string& pssmName, const vector<string> & PSSMLines, const vector<char> & correspondingCharacters, const vector<double> & aaFreq, size_t pseudoCountValue) :
	_Nsite(numberOfSupportingSequences), _PSSM_name(pssmName), _correspondingCharacters(correspondingCharacters), _aaFreq(aaFreq) {
	size_t numberOfPositions = PSSMLines.size(); // init lines of matrix
	_PSSMmatrix.resize(numberOfPositions);
	for (size_t MotifPos = 0; MotifPos < PSSMLines.size(); MotifPos++) // fill in lines
	{
		string line = PSSMLines[MotifPos]; //s = "0 0.926111458985598 0.000626174076393237 0 0.00438321853475266 0.0050093926111459 0 0.00125234815278647 0.00187852222917971 0.00250469630557295 0.00250469630557295 0 0.000626174076393237 0 0.00438321853475266 0.0219160926737633 0.00751408891671885 0 0.0137758296806512 0.00751408891671885";
		stringstream line_element(line); // split line by spaces
		istream_iterator<string> begin(line_element);
		istream_iterator<string> end;
		vector<string> vstrings(begin, end);
		//		_PSSMmatrix[MotifPos].resize(vstrings.size()+1); // last position is for padding gap character
		_PSSMmatrix[MotifPos].resize(vstrings.size());
		for (size_t i = 0; i < vstrings.size(); i++)
		{
			double lol = atof(vstrings[i].c_str());
			_PSSMmatrix[MotifPos][i] = lol;
		}
		//		_PSSMmatrix[MotifPos][vstrings.size()]=0; // last position is for padding gap character
	}
	_PSSM_consensus_seq = computeConsensusSeq();
	//	cout<<"done with pssm matrix before pseudocounts"<<endl;
	//	cout<<"num elements in first row for example is "<<_PSSMmatrix[0].size()<<endl;
	//	cout<<"first row is:"<<endl;
	double sum_first_row = 0;
	for (size_t i = 0; i < _PSSMmatrix[0].size(); i++)
	{
		sum_first_row = sum_first_row + _PSSMmatrix[0][i];
		//		cout<<_PSSMmatrix[0][i]<<endl;
	}
	//	cout<<"sum of elements in first row for example is "<<sum_first_row<<endl;
	//	cout<<"gap element in first row for example is "<<_PSSMmatrix[0][_PSSMmatrix[0].size()-1]<<endl;

	addPseudoCount(pseudoCountValue);
	for (size_t i = 0; i < _PSSMmatrix.size(); ++i) {
		vector<double> tmpLogValues;
		for (size_t k = 0; k < alphabetSize(); ++k) {
			tmpLogValues.push_back(log(_PSSMmatrix[i][k]));
		}
		_PSSMmatrixLOG.push_back(tmpLogValues);
	}

	_cummulative_aaFreq.push_back(_aaFreq[0]);
	for (size_t j = 1; j < _aaFreq.size(); ++j) {
		_cummulative_aaFreq.push_back(_aaFreq[j] + _cummulative_aaFreq[j - 1]);
	}
}



void PSSM::addPseudoCount(size_t val)
{
	if (_Nsite == 0)
	{
		cerr << "ERROR: Cna't employ AddPseudoCount on motif: '" << _PSSM_name << "' because Nsite is not defiend..." << endl;
		return;
	}
	for (size_t MotifPos = 0; MotifPos < _PSSMmatrix.size(); MotifPos++) {
		for (size_t i = 0; i < _PSSMmatrix[MotifPos].size(); i++)
		{
			double currValue = _PSSMmatrix[MotifPos][i];
			double newVal = ((currValue*_Nsite) + val) / (_Nsite + (val*_PSSMmatrix[MotifPos].size()));
			//			newVal=log(newVal)/log(2.0);
			_PSSMmatrix[MotifPos][i] = newVal;
		}
	}
	//	cout<<"num seqs was: "<<_Nsite<<endl;
	//	cout<<"done with pssm matrix after pseudocounts"<<endl;
	//	cout<<"num elements in first row for example is "<<PSSMmatrix[0].size()<<endl;
	//	cout<<"first row is:"<<endl;
	//	double sum_first_row = 0;
	//	for(int i = 0; i < PSSMmatrix[0].size(); i++)
	//	{
	//		sum_first_row = sum_first_row + PSSMmatrix[0][i];
	//		cout<<PSSMmatrix[0][i]<<endl;
	//	}
	//	cout<<"sum of elements in first row for example is "<<sum_first_row<<endl;
	//	cout<<"gap element in first row for example is "<<PSSMmatrix[0][PSSMmatrix[0].size()-1]<<endl;
}


double getVectorPercentile(const vector<double>& vec, double quantile) {
	vector<double> sortedVec = vec;
	sort(sortedVec.begin(), sortedVec.end());
	size_t placeValue = static_cast<size_t>(vec.size()*quantile);
	return sortedVec[placeValue];
}


double PSSM::getScoreOf_a_RandomPeptide(const vector<size_t> & randomPeptideEntries,size_t peptideSize) {
	
	double maxScore = -10000000.0;
	for (size_t i = 0; i < _PSSMmatrix.size() + peptideSize - 1; ++i) {
		double totalScore = 0.0;
		if (i < peptideSize) {
			for (size_t k = peptideSize-i-1; k < peptideSize; ++k) {//  // k here is the position in the peptide
				totalScore += getLogValue(k, randomPeptideEntries[k]);// in this specific case only, k is also the position in the PSSM...
			}
			totalScore += ((peptideSize - i - 1)*log(1.0 / alphabetSize())); // giving a penalty of 1/alphabetsize for each position that is not aligned to the PSSM
		}
		else if (i >= _PSSMmatrix.size()) {
			for (size_t k =0; k < peptideSize-i+_PSSMmatrix.size(); ++k) {  // k here is the position in the peptide
				totalScore += getLogValue(i-peptideSize+k, randomPeptideEntries[k]);
			}
			totalScore += ((i - _PSSMmatrix.size() + 1)*log(1.0 / alphabetSize()));
		}
		else {// peptide is within the PSSM
			for (size_t k = 0; k <peptideSize; ++k) { // k here is the position in the peptide
				totalScore += getLogValue(i- peptideSize+k+1, randomPeptideEntries[k]);
			}
		}
		if (totalScore > maxScore) maxScore = totalScore;
	}
	return maxScore;
}



double PSSM::giveScorePercentile(size_t peptideSize, double percentile, size_t numberOfSimulations) {
	mt19937 gen(0); //Standard mersenne_twister_engine seeded with rd()
	uniform_real_distribution<> dis(0, 1);

	vector<double> simulationsResults;
	for (size_t i = 0; i < numberOfSimulations; ++i) {
		
		// generating the random number vector;
		vector<size_t> randomPeptideEntries;
		for (size_t j = 0; j < peptideSize; ++j) {
			double randNumBetweenZeroAndOne = dis(gen);
			for (size_t k = 0; k < _cummulative_aaFreq.size(); ++k) {
				if (randNumBetweenZeroAndOne < _cummulative_aaFreq[k]) {
					randomPeptideEntries.push_back(k);
					break;
				}
			}
		}
		
		double aSingleSimulationResult = getScoreOf_a_RandomPeptide(randomPeptideEntries, peptideSize);
		simulationsResults.push_back(aSingleSimulationResult);
	}
	double resPercentile = getVectorPercentile(simulationsResults, percentile);
	return resPercentile;
}