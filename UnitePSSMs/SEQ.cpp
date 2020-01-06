#include "SEQ.h"
 /// NOTE: to suppurt regex, gcc>=4.9 is needed!!!
void SEQ::setCopyNmber()
{
		_CopyNumber=1.0;
		smatch match;
		regex rgx(".*Repeats_([\\d.]+)_Type");
		if (regex_search(_Seq_Name, match, rgx))
		{
			string CopyNumber_str=match[1].str();
			_CopyNumber=atof(CopyNumber_str.c_str());
		}
}
void SEQ::setSeq_Type()
{
	int seq_length=_Seq.length();
	stringstream ss;
	if (seq_length<1)
	{
		cout<<"Invalid sequence "<<_Seq_Name<<" ignored..."<<endl;
	}
	char first = _Seq[0];
	char last = _Seq[seq_length-1];
	if ((first == 'C') && (last == 'C')) 
	{
		seq_length=seq_length-2;
		ss<<"C"<<seq_length<<"C";
	}
	else
	{
		ss<<seq_length;
	}
	_Seq_Type=ss.str();
}
