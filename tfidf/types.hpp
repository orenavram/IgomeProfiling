#pragma once
#include <string>
#include <map>
#include <vector>
using namespace std;

// types
enum TFMethod {
    TFMethod_NONE,
    TFMethod_BOOL,
    TFMethod_TERMS,
    TFMethod_LOG,
    TFMethod_AUGMENT
};

class Meme;
class MemeSample;
typedef map<string, int> SequencesCount;
typedef map<string, string> SamplesLabel;
typedef map<string, MemeSample*> MemeSampleMap;
typedef map<string, Meme*> MemesMap;
typedef map<string, string> SamplesBC;
typedef vector<string> MemesList;
