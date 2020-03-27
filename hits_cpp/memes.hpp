#pragma once
#include "types.hpp"

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
