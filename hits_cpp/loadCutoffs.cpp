#include <iostream>
#include <string>
#include <fstream>
#include <regex>

#include "types.hpp"
#include "meme.hpp"
#include "memes.hpp"

using namespace std;
void loadCutoffs(string cutoffsPath, Memes& memes, int limit = 0, bool verbose = false) {
    ifstream file(cutoffsPath);
    string line;

    regex pattern("([^\\t]+)\\t([^\\t]+)");
    smatch matches;
    Meme* meme;
    int total = 0;

    while (getline(file, line)) {
        regex_search(line, matches, pattern);
        if (matches[1].compare("###") == 0) {
            if (limit && total == limit) {
                break;
            }
            meme = &memes.getMemes()[matches[2]];
            total++;
        } else {
            if (verbose) {
                cout << matches[1] << ": " << matches[2] << endl;
            }
            meme->getCuttofs()[matches[1]] = stod(matches[2]);
        }
    }
}