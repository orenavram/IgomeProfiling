#include <string>
#include <iostream>
#include <fstream>

#include "types.hpp"
#include "scans.hpp"
#include "meme.hpp"
#include "memeSample.hpp"

using namespace std;
void writeResults(MemesList& memes, Scans& scans, string& outputPath, string& bc) {
    ofstream hitsFile(outputPath + bc + "_hits.csv");
    ofstream scoreFile(outputPath + bc + "_values.csv");
    auto header = "sample_name,label";
    hitsFile << header;
    scoreFile << header;

    // Transpose for easier write
    auto samplesCount = scans.getSamplesLabel().size();
    auto memesCount = scans.getMemes().size();
    auto hitsValues = new double[memesCount * samplesCount];
    auto scoreValues = new double[memesCount * samplesCount];
    
    auto memesIter = memes.begin();
    auto memesEnd = memes.end();
    int memeIndex = 0;
    while (memesIter != memesEnd) {
        hitsFile << "," << *memesIter;
        scoreFile << "," << *memesIter;

        auto meme = scans.getMemes().find(*memesIter);

        auto samplesIter = meme->second->getSamples().begin();
        auto samplesEnd = meme->second->getSamples().end();
        int sampleIndex = 0;
        while (samplesIter != samplesEnd) {
            auto hits = samplesIter->second->getHitsCount();
            auto score = samplesIter->second->getScore();
            hitsValues[sampleIndex * memesCount + memeIndex] = hits;
            scoreValues[sampleIndex * memesCount + memeIndex] = score;
            sampleIndex++;
            samplesIter++;
        }
        memeIndex++;
        memesIter++;
    }
    hitsFile << endl;
    scoreFile << endl;

    auto samplesIter = scans.getSamplesLabel().begin();
    auto samplesEnd = scans.getSamplesLabel().end();
    int sampleIndex = 0;
    while (samplesIter != samplesEnd) {
        hitsFile << samplesIter->first << "," << samplesIter->second;
        scoreFile << samplesIter->first << "," << samplesIter->second;
        
        auto memesIter = memes.begin();
        auto memesEnd = memes.end();
        int memeIndex = 0;
        while (memesIter != memesEnd) {
            hitsFile << "," << hitsValues[sampleIndex * memesCount + memeIndex];
            scoreFile << "," << scoreValues[sampleIndex * memesCount + memeIndex];
            memeIndex++;
            memesIter++;
        }
        hitsFile << endl;
        scoreFile << endl;
        sampleIndex++;
        samplesIter++;
    }
}
