#include <string>
#include <iostream>
#include <random>

#include "shufflePatterns.hpp"
#include "meme.hpp"

#include "memes.hpp"

using namespace std;

Memes loadMemes(string memePath, int limit, bool verbose);

class MemeShuffler {
public:
    MemeShuffler(Meme& source, int maxShuffles, ShufflePatterns& patterns) :
        _source(source), _patterns(patterns), _currentShuffle(-1) {
        auto len = patterns.size();
        if (len < maxShuffles) {
            maxShuffles = len;
        }
        this->_maxShuffles = maxShuffles;
    }

    bool next() {
        _currentShuffle++;
        return hasNext();
    }

    Meme generate() {
        auto meme = Meme(this->_source);
        meme.getRows().clear();

        auto pattern = this->_patterns[this->_currentShuffle];
        auto iter = pattern.begin();
        auto end = pattern.end();
        while (iter != end) {
            meme.getRows().push_back(this->_source.getRows()[*iter]);
            iter++;
        }
        return meme;
    }

    bool hasNext() {
        return _currentShuffle < _maxShuffles;
    }

    int getMaxShuffles() {
        return this->_maxShuffles;
    }
private:
    Meme _source; //TODO pointer?
    ShufflePatterns _patterns; //TODO pointer?
    int _maxShuffles;
    int _currentShuffle;
};


vector<int>* getPattern(int length) {
    auto iter = _shuffle_sequences.find(length);
    if (iter == _shuffle_sequences.end()) {
        return nullptr;
    }
    auto pattern = &iter->second[1];
    return pattern;
}

MemeShuffler getShuffler(Meme& meme, int maxPatterns) {
    auto length = meme.getRows().size();
    auto iter = _shuffle_sequences.find(length);
    if (iter == _shuffle_sequences.end()) {
        ShufflePatterns patterns;
        return MemeShuffler(meme, -1, patterns);
    }
    return MemeShuffler(meme, maxPatterns, iter->second);
}

// int main() {
//     // cout << "hello" << endl;
//     // auto pattern = getPattern(7);
//     // cout << (*pattern)[0] << endl;

//     Memes memes = loadMemes("../output/mock_small_pval/analysis/motif_inference/17b/meme.txt", 0, false);
//     auto memesIter = memes.getMemes().begin();
//     cout << memesIter->first << endl;
//     auto meme = memesIter->second;
//     cout << "rows: " << meme.getRows().size() << ", a: " << meme.getALength() << endl;
//     auto shuffler = getShuffler(meme, 5);
//     while (shuffler.next()) {
//         auto newMeme = shuffler.generate();
//         cout << newMeme.getMotif() << endl;
//         cout << newMeme.getRows().size() << endl;
//     }
// }