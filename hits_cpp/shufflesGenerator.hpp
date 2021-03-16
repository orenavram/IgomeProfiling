#pragma once

#include "meme.hpp"

using namespace std;

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
        //check that we don't change place of the same data row:
        for (int i = 0; i < this->_source.getRows().size(); i++){
            if (this->_source.getRows()[i] == meme.getRows()[i]){
                if (this->_source.getRows().size() > _maxShuffles){
                    _currentShuffle++;
                    _maxShuffles++;
                    auto meme_new = this->generate();
                    return meme_new;
                }    
            }
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

MemeShuffler getShuffler(Meme& meme, int maxPatterns);
