#ifndef SEQ_____H
#define SEQ_____H

#include <vector>
#include <string>
#include <regex>
#include <iostream>
#include <sstream>      // std::stringstream

using namespace std;


class SEQ {
public:
	string _Seq;
	string _Seq_Name;
	string _Seq_Type;
	double _CopyNumber;
//	double _HitScore;
//	SEQ(){_CopyNumber=0;setSeq("");setHitScore(-1);};
	SEQ(){_CopyNumber=0;setSeq("");};
	SEQ(const string & SeqString,const string & SeqName, const double CopyNumber): _Seq(SeqString), _Seq_Name(SeqName), _CopyNumber(CopyNumber)
	{
	};
	~SEQ(){};
	void setName(const string name ){ _Seq_Name = name;}
	void setCopyNumber (const double copies){_CopyNumber=copies;}
	void setSeq(const string seq){_Seq=seq;}
//	void setHitScore (const double score){_HitScore=score;} // TO DO make it nicer...
	void setCopyNmber ();
	void setSeq_Type();

};

#endif