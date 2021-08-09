#ifndef SEQ_____H
#define SEQ_____H

#include <vector>
#include <string>
#include <regex>
#include <iostream>
#include <sstream>      // std::stringstream

using namespace std;
#include "alphabet.h"

class SEQ {
public:
	vector<size_t> _Seq;
	string _Seq_Name;
	string _Seq_Type;
	double _CopyNumber;
	//SEQ(){_CopyNumber=0;};
	SEQ(const vector<size_t> & SeqString,const string & SeqName, const double CopyNumber, alphabet& alph): _Seq(SeqString), _Seq_Name(SeqName), _CopyNumber(CopyNumber), _alph(alph) {
		setSeq_Type();
		setCopyNmber();
	};
	SEQ(string & SeqString, const string & SeqName, const double CopyNumber, alphabet& alph, bool isSetCopyNumber);

	~SEQ(){};
	void setName(const string name ){ _Seq_Name = name;}
	void setCopyNumber (const double copies){_CopyNumber=copies;}
	void setSeq(const vector<size_t>& seq){_Seq=seq;}
//	void setHitScore (const double score){_HitScore=score;} // TO DO make it nicer...
	void setCopyNmber ();
	void setSeq_Type();
	alphabet& _alph;
	string getStringOfSeq() const;
	string getStringOfSeq(size_t i) const; // print from position i till the end
};

#endif