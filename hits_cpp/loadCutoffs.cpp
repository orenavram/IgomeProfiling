#include <iostream>
#include <string>
#include <fstream>
#include <regex>

#include "types.hpp"
#include "meme.hpp"
#include "memes.hpp"

using namespace std;
void loadCutoffs(string cutoffsPath, Memes& memes) {
    ifstream file(cutoffsPath);
    string line;

    regex pattern("([^\\t]+)\\t([^\\t]+)");
    smatch matches;
    Meme* meme;

    while (getline(file, line)) {
        regex_search(line, matches, pattern);
        if (matches[1].compare("###") == 0) {
            meme = &memes.getMemes()[matches[2]];
        } else {
            cout << matches[2] << endl;
            meme->getCuttofs()[matches[1]] = stod(matches[2]);
        }
    }
}