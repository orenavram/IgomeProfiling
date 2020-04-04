#include <math.h>
#include "types.hpp"
#include "scans.hpp"
#include "meme.hpp"
#include "memeSample.hpp"

void calculateScores(Scans& scans, TFMethod eMethod, float augmentFactor, bool isUseAllSequences) {
    auto memes = &scans.getMemes();
    auto sequencesCount = &scans.getSequencesCount();

    auto memesIter = memes->begin();
    auto memesEnd = memes->end();
    auto bcSequencesEnd = scans.getBCSequences().end();
    double memesCount = (double)memes->size();
    // cout << "memes count: " << memesCount << endl;
    while (memesIter != memesEnd) {
        auto samplesIter = memesIter->second->getSamples().begin();
        auto samplesEnd = memesIter->second->getSamples().end();
        while (samplesIter != samplesEnd) {
            auto sequencesIter = samplesIter->second->getSequences().begin();
            auto sequencesEnd = samplesIter->second->getSequences().end();
            double score = 0;
            while (sequencesIter != sequencesEnd) {
                if (isUseAllSequences || scans.getBCSequences().find(sequencesIter->first) != bcSequencesEnd) {
                    // TODO calculate cosine similarity with all relevant mems if on
                    auto sequenceMemesCount = sequencesCount->find(sequencesIter->first)->second;
                    double idf = log(memesCount / sequenceMemesCount);
                    double tf = 0;
                    switch (eMethod) {
                    case TFMethod_BOOL:
                        tf = 1;
                        break;
                    case TFMethod_TERMS:
                        tf = (double)sequencesIter->second / samplesIter->second->getHitsCount();
                        break;
                    case TFMethod_LOG:
                        tf = log(1 + sequencesIter->second);
                        break;
                    case TFMethod_AUGMENT:
                        tf = augmentFactor + (1 - augmentFactor) * 
                            ((double)sequencesIter->second / samplesIter->second->getMaxSequenceCount());
                        break;
                    }
                    score += tf * idf;
                }
                sequencesIter++;
            }
            samplesIter->second->setScore(score);
            samplesIter++;
        }
        memesIter++;
    }
}
