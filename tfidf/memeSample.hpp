#pragma once
#include "types.hpp"

class MemeSample {
public:
    MemeSample(int hitscount) :
        _hitsCount(hitscount),
        _score(0),
        _maxSequenceCount(0) {
    }

    void addSequence(string& sequence, int count) {
        this->_sequences[sequence] = count;
        if (count > this->_maxSequenceCount) {
            this->_maxSequenceCount = count;
        }
    }

    int getHitsCount() {
        return this->_hitsCount;
    }

    SequencesCount& getSequences() {
        return this->_sequences;
    }

    int getMaxSequenceCount() {
        return this->_maxSequenceCount;
    }

    double getScore() {
        return this->_score;
    }

    void setScore(double score) {
        this->_score = score;
    }
private:
    int _hitsCount;
    SequencesCount _sequences;
    int _maxSequenceCount;
    double _score;
};
