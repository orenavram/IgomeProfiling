#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <math.h>
#include <fstream>
#include <regex>
#include <chrono>

#include "cxxopts.hpp"

// Alphabet of 20 + gap
#define MAX_MEME_COLUMNS 21

using namespace std;
class Meme;
typedef map<string, Meme> MemesMap;
typedef map<string, list<string>*> SequencesMap;
typedef map<char, int> AlphabetMap;
typedef map<string, double> CutoffsMap;
typedef vector<vector<double>> MemeRows;

string& ltrim(string& str, const string& chars = "\t\n\v\f\r ")
{
    str.erase(0, str.find_first_not_of(chars));
    return str;
}
 
string& rtrim(string& str, const string& chars = "\t\n\v\f\r ")
{
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}
 
string& trim(string& str, const string& chars = "\t\n\v\f\r ")
{
    return ltrim(rtrim(str, chars), chars);
}

class Meme {
public:
    Meme() {

    }

    string getMotif() {
        return this->_motif;
    }

    void setMotif(string motif) {
        this->_motif = motif;
    }

    int getALength() {
        return this->_alength;
    }

    void setALength(int alength) {
        this->_alength = alength;
    }

    int getNSites() {
        return this->_nsites;
    }

    void setNSites(int nsites) {
        this->_nsites = nsites;
    }

    MemeRows& getRows() {
        return this->_rows;
    }

    CutoffsMap& getCuttofs() {
        return this->_cutoffs;
    }

    void normalize() {
        if (this->_nsites == 0) {
            return;
        }
        auto rows = this->_rows.size();
        auto columns = this->_alength + 1;
        
        auto iter = this->_rows.begin();
        auto end = this->_rows.end();

        for (auto rowsIter = this->_rows.begin(); rowsIter != this->_rows.end(); ++rowsIter) {
            for (auto columnsIter = rowsIter->begin(); columnsIter != rowsIter->end(); ++columnsIter) {
                *columnsIter = ((*columnsIter * this->_nsites) + 1) / (this->_nsites + columns);
            }
        }
    }
private:
    string _motif;
    int _alength;
    int _nsites;
    MemeRows _rows;
    CutoffsMap _cutoffs;
};

class Memes {
public:
    Memes(AlphabetMap alphabet, MemesMap memes) :
        _alphabet(alphabet), _memes(memes) {
    }

    AlphabetMap& getAlphabet() {
        return this->_alphabet;
    }

    MemesMap& getMemes() {
        return this->_memes;
    }
private:
    AlphabetMap _alphabet;
    MemesMap _memes;
};

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

void loadCutoffs(string cutoffsPath, Memes& memes) {
    ifstream file(cutoffsPath);
    string line;

    regex pattern("([^\\t]+)\\t([^\\t]+)");
    smatch matches;
    Meme* meme;

    while (getline(file, line)) {
        regex_search(line, matches, pattern);
        if (matches[1].compare("###") == 0) {
            meme = &memes.getMemes()[matches[2]]; // TODO does this copy?
        } else {
            cout << matches[2] << endl;
            meme->getCuttofs()[matches[1]] = stod(matches[2]);
        }
    }
}

SequencesMap loadSequences(string faaPath) {
    SequencesMap sequences;
    ifstream file(faaPath);
    string line;
    list<string>* sequencesByType;
    
    regex pattern("^>.+_(.+)$");
    smatch matches;
    auto end = sequences.end();

    int count = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
            auto lastIndex = line.find_last_of('_');
            auto seqType = line.substr(lastIndex + 1);
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
            count++;
        }
    }

    cout << "total sequences: " << count << endl;
    return sequences;
}

bool isHit(Meme& meme, AlphabetMap& alphabet, string seqType, string& seq) {
    auto iter = meme.getCuttofs().find(seqType);
    if (iter == meme.getCuttofs().end()) {
        return false;
    }
    auto cutoffValue = iter->second;
    auto seqLen = seq.length();
    auto memeLen = meme.getRows().size(); // TODO get by variable?

    int start = -memeLen + 1;
    int end = min(memeLen, seqLen);
    double totalScore = 0;
    double score = 0;
    char c;
    for (int i = start; i < end; i++) {
        totalScore = 1;
        for (int j = 0; j < memeLen; j++) {
            c = ((i + j < 0) || (i + j) >= seqLen) ? '-' : seq[i + j];
            score = meme.getRows()[j][alphabet[c]];
            // TODO check if equal to -inf, set total_score and break
            totalScore *= score;
        }
        if (log(totalScore) > cutoffValue) {
            return true;
        }
    }
    return false;
}

int getHits(Memes& memes, SequencesMap sequences) {
    cout << "GET HITS" << endl;
    auto alphabet = memes.getAlphabet();
    auto memesIter = memes.getMemes().begin();
    auto memesEnd = memes.getMemes().end();

    int hits = 0;
    int counter = 0;
    while (memesIter != memesEnd) {
        auto sequencesTypesIter = sequences.begin();
        auto sequencesTypesEnd = sequences.end();
        while (sequencesTypesIter != sequencesTypesEnd) {
            auto sequencesIter = sequencesTypesIter->second->begin();
            auto sequencesEnd = sequencesTypesIter->second->end();
            while (sequencesIter != sequencesEnd) {
                if ((counter % 10000) == 0) {
                    cout << "seq: " << counter << ", hits: " << hits << endl;
                }
                counter++;
                if (isHit(memesIter->second, alphabet, sequencesTypesIter->first, *sequencesIter)) {
                    hits++;
                }
                sequencesIter++;
            }
            sequencesTypesIter++;
        }
        memesIter++;
    }
    cout << "total hits: " << hits << endl;
    return hits;
}

int main(int argc, char *argv[])
{
    // cxxopts::Options options("MyProgram", "One line description of MyProgram");
    // options.add_options()("d,debug", "Enable debugging")("f,file", "File name", cxxopts::value<string>());
    // auto result = options.parse(argc, argv);
    // cout << result["file"].as<string>() << endl;

    string memesPath = "/home/shacharm/test2/01_exp4+9_exp4_naive_meme_20.txt";
    string cutoffsPath = "/home/shacharm/test2/02_exp4+9_exp4_naive_cutoff_20.txt";
    string sequencesPath = "/home/shacharm/test2/03_exp4+9_exp4_naive_c.faa";

    auto begin = chrono::steady_clock::now();

    Memes memes = loadMemes(memesPath);
    loadCutoffs(cutoffsPath, memes);
    SequencesMap sequences = loadSequences(sequencesPath);
    getHits(memes, sequences);

    auto end = chrono::steady_clock::now();
    cout << chrono::duration_cast<chrono::seconds>(end - begin).count() << "[s]" << endl;
}
