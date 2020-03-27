#include "utils.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <regex>

using namespace std;
void tolower(string& str) {
    for_each(str.begin(), str.end(), [](char & c){
	    c = ::tolower(c);
    });
}

TFMethod parseMethod(string method) {
     tolower(method);

     if (method.compare("boolean") == 0) {
         return TFMethod_BOOL;
     }

     if (method.compare("terms") == 0) {
         return TFMethod_TERMS;
     }

     if (method.compare("log") == 0) {
         return TFMethod_LOG;
     }

     if (method.compare("augmented") == 0) {
         return TFMethod_AUGMENT;
     }

     return TFMethod_NONE;
}

SamplesBC loadSamplesToBC(string samples2bcPath) {
    SamplesBC samplesToBC;
    ifstream file(samples2bcPath);

    string line;
    smatch matches;
    regex pattern("^([^#\t]+?)\\t(.+?)$");
    while (getline(file, line)) {
        if (regex_search(line, matches, pattern)) {
            // cout << matches[1] << ": " << matches[2] << endl;
            samplesToBC[matches[1]] = matches[2];
        }
    }
    return samplesToBC;
}
