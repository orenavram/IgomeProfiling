#include <string>
#include <iostream>
#include <fstream>
#include <regex>

#include "types.hpp"
#include "meme.hpp"
#include "memes.hpp"
#include "trim.hpp"

using namespace std;

Memes loadMemes(string memePath) {
    ifstream file(memePath);
    string line;

    AlphabetMap alphabet;
    MemesMap memes;

    regex motifPattern("alength= (\\d+) w= (\\d+) nsites= (\\d+)");
    regex rowPattern("(\\d+\\.?\\d*)");
    auto regexEnd = sregex_iterator();
    smatch matches;
    int alength;
    int rows;

    while (getline(file, line)) {
        if (line.rfind("MOTIF", 0) == 0) {
            Meme meme;
            string motif = line.substr(6);
            meme.setMotif(rtrim(motif));
            
            getline(file, line);
            regex_search(line, matches, motifPattern);
            alength = stoi(matches[1]);
            meme.setALength(alength);
            rows = stoi(matches[2]);
            meme.setNSites(stoi(matches[3]));

            for (int i = 0; i < rows; i++) {
                getline(file, line);
                vector<double> row;
                sregex_iterator iter(line.begin(), line.end(), rowPattern);
                while (iter != regexEnd) {
                    row.push_back(stod(iter->str()));
                    iter++;
                }
                row.push_back(0); // gap
                meme.getRows().push_back(row);
            }
            meme.normalize();
            memes[motif] = meme;
        } else if (line.rfind("ALPHABET=", 0) == 0) {
            int start = 10;
            int end = line.length() - start;
            for (int i = 0; i < end; i++) {
                alphabet[line[start + i]] = i;
            }
            alphabet['-'] = end;
        }
    }
    return Memes(alphabet, memes);
}
