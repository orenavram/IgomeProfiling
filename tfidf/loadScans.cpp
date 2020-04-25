#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <stdio.h>
#include <dirent.h>

#include "types.hpp"
#include "scans.hpp"
#include "meme.hpp"
#include "memeSample.hpp"
#include "utils.hpp"

using namespace std;
void loadScan(Scans& scans, string& scanPath, SamplesBC& samplesToBC, string& bc, 
    bool isMulticlass, bool readRankAsScore) {
    ifstream file(scanPath);

    auto memes = &scans.getMemes();
    auto samplesLabel = &scans.getSamplesLabel();
    auto sequencesCount = &scans.getSequencesCount();

    regex filePattern("^.+/(.+?)_peptides_vs_(.+?)_motifs_");
    regex motifPattern("^MOTIF (.+?)_");
    regex hitsPattern("^([^ ]+) (.+)$");

    smatch matches;
    regex_search(scanPath, matches, filePattern);
    if (bc.compare(matches[2]) != 0) {
        return; // scan is not part of biological condition
    }

    string other = "other";
    string sample = matches[1];
    auto isSampleOfBC = bc.compare(samplesToBC[sample]) == 0;

    if (isMulticlass) {
        samplesLabel->insert(pair<string, string>(sample, samplesToBC[sample]));
    } else if (isSampleOfBC) {
        samplesLabel->insert(pair<string, string>(sample, bc));
    } else {
        samplesLabel->insert(pair<string, string>(matches[1], other));
    }

    string line;
    Meme* meme;
    MemeSample* memeSample;
    while (getline(file, line)) {
        if (regex_search(line, matches, motifPattern)) {
            string motif = matches[1];
            auto memeIter = memes->find(motif);
            if (memeIter == memes->end()) {
                meme = new Meme();
                memes->insert(pair<string, Meme*>(motif, meme));
            } else {
                meme = memeIter->second;
            }
            
            getline(file, line); // hits
            regex_search(line, matches, hitsPattern);
            int hits = stoi(matches[2]);
            memeSample = new MemeSample(hits);
            // cout << "motif " << motif << ", sample hits: " << sample << " = " << hits << endl;
            if (readRankAsScore) {
                getline(file, line); // shuffles
                getline(file, line); // rank
                regex_search(line, matches, hitsPattern);
                auto rank = stod(matches[2]);
                memeSample->setScore(rank);
            }
            meme->addSample(sample, memeSample);
        }
        else if (regex_search(line, matches, hitsPattern)) // sequences
        {
            string sequence = matches[1];
            int hits = stoi(matches[2]);
            
            memeSample->addSequence(sequence, hits);

            auto sequencesCountIter = sequencesCount->find(sequence);
            if (sequencesCountIter == sequencesCount->end()) {
                sequencesCount->insert(pair<string, int>(sequence, 1));
            } else {
                sequencesCountIter->second += 1;
            }

            if (isSampleOfBC) {
                scans.getBCSequences().insert(sequence);
            }
        }
    }
}

Scans loadScans(string scanPath, string samples2bcPath, string& bc, 
    bool isMultiClass, bool readRankAsScore) {
    Scans scans;
    auto samplesToBC = loadSamplesToBC(samples2bcPath);

    struct dirent *entry = nullptr;
    DIR *dp = nullptr;
    dp = opendir(scanPath.c_str());
    if (dp != nullptr) {
        while ((entry = readdir(dp))) {
            string entryName(entry->d_name);
            string entryPath = scanPath + entryName;
            // cout << entryPath << endl;
            loadScan(scans, entryPath, samplesToBC, bc, isMultiClass, readRankAsScore);
        }
    }
    closedir(dp);
    
    return scans;
}
