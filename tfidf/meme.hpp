#pragma once
#include "types.hpp"

class Meme {
public:
    Meme() {
    }

    void addSample(string& sample, MemeSample* memeSample) {
        this->_samples[sample] = memeSample;
    }

    MemeSampleMap& getSamples() {
        return this->_samples;
    }
private:
    MemeSampleMap _samples;
};
