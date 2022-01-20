#include <iostream>
#include <string>
#include <fstream>
#include <regex>
#include <sstream>
#include <vector>

#include "types.hpp"
#include "trim.hpp"

using namespace std;
SequencesMap loadSequences(string faaPath, int& numSequences, SequencesRpmMap& sequncesRpm, bool useRpmFaaScanning, bool verbose) {
    SequencesMap sequences;
    ifstream file(faaPath);
    string line;
    list<string>* sequencesByType;

    regex pattern("^>.+_(.+)$");
    smatch matches;
    auto end = sequences.end();
    int count = 0;
    int uniquePeptides = 1;
    string seqType;
    while (getline(file, line)) {
        if (line[0] == '>') {
            if (useRpmFaaScanning) {
                replace(line.begin(), line.end(), '_', ' ');  // replace '_' by ' '
                istringstream iss(line);
                vector<string> result;
                for (string s;iss>>s;) {
                   result.push_back(s);
                }
                seqType = result[3];
                uniquePeptides = stoi(result[9]);
            } else {
                auto lastIndex = line.find_last_of('_');
                seqType = line.substr(lastIndex + 1);
            }
            rtrim(seqType);
            auto iter = sequences.find(seqType);
            if (iter == end) {
                sequencesByType = new list<string>();
                sequences[seqType] = sequencesByType;
            } else {
                sequencesByType = iter->second;
            }
            
            getline(file, line);
            sequencesByType->push_back(line);
            if (useRpmFaaScanning) {
                sequncesRpm[line] = uniquePeptides;
                count += uniquePeptides;
            } else {
                auto iter = sequncesRpm.find(line);
                if (iter == sequncesRpm.end()) {
                    sequncesRpm[line] = 1;
                } else {
                    iter->second++;
                }
                count++;
            }
        }
    }
    numSequences = count;
    cout << "total sequences: " << count << endl;
    return sequences;
}
