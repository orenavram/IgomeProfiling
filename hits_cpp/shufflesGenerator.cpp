#include <string>
#include <iostream>

#include "shufflesGenerator.hpp"
#include "shufflePatterns.hpp"

#include "memes.hpp"

using namespace std;

Memes loadMemes(string memePath, int limit, bool verbose);

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