#pragma once
#include "types.hpp"

class Scans {
public:
    Scans() {
    }

    MemesMap& getMemes() {
        return this->_memes;
    }

    SamplesLabel& getSamplesLabel() {
        return this->_samplesLabel;
    }

    SequencesCount& getSequencesCount() {
        return this->_sequencesMemeCount;
    }
private:
    MemesMap _memes;
    SamplesLabel _samplesLabel;
    SequencesCount _sequencesMemeCount;
};
