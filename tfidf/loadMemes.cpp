#include <string>
#include <iostream>
#include <fstream>
#include <regex>

#include "types.hpp"

using namespace std;

MemesList loadMemes(string memesPath) {
    MemesList list;
    
    ifstream file(memesPath);
    string line;

    regex pattern("MOTIF (.+?)_");
    smatch matches;

    while (getline(file, line)) {
        if (regex_search(line, matches, pattern)) {
            list.push_back(matches[1]);
        }
    }

    return list;   
}