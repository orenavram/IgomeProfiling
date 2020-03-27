#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <array>
#include <regex>
#include <stdio.h>
#include <dirent.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>

#include "cxxopts.hpp"

using namespace std;

// types
enum TFMethod {
    TFMethod_NONE,
    TFMethod_BOOL,
    TFMethod_TERMS,
    TFMethod_LOG,
    TFMethod_AUGMENT
};

class Meme;
class MemeSample;
typedef map<string, int> SequencesCount;
typedef map<string, string> SamplesLabel;
typedef map<string, MemeSample*> MemeSampleMap;
typedef map<string, Meme*> MemesMap;
typedef map<string, string> SamplesBC;

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

// Utils
void tolower(string& str) {
    for_each(str.begin(), str.end(), [](char & c){
	    c = ::tolower(c);
    });
}

TFMethod parseMethod(string method) {
     tolower(method);

     if (method.compare("boolean") == 0) {
         return TFMethod_BOOL;
     }

     if (method.compare("terms") == 0) {
         return TFMethod_TERMS;
     }

     if (method.compare("log") == 0) {
         return TFMethod_LOG;
     }

     if (method.compare("augmented") == 0) {
         return TFMethod_AUGMENT;
     }

     return TFMethod_NONE;
}

// Methods
SamplesBC loadSamplesToBC(string samples2bcPath) {
    SamplesBC samplesToBC;
    ifstream file(samples2bcPath);

    string line;
    smatch matches;
    regex pattern("^([^#\t]+?)\\t(.+?)$");
    while (getline(file, line)) {
        if (regex_search(line, matches, pattern)) {
            // cout << matches[1] << ": " << matches[2] << endl;
            samplesToBC[matches[1]] = matches[2];
        }
    }
    return samplesToBC;
}

void loadScan(Scans& scans, string& scanPath, SamplesBC& samplesToBC, string& bc, bool isMulticlass) {
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
    if (isMulticlass) {
        samplesLabel->insert(pair<string, string>(sample, samplesToBC[sample]));
    } else if (bc.compare(samplesToBC[sample]) == 0) {
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
        }
    }
}

Scans loadScans(string scanPath, string samples2bcPath, string& bc, bool isMultiClass) {
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
            loadScan(scans, entryPath, samplesToBC, bc, isMultiClass);
        }
    }
    closedir(dp);
    
    return scans;
}

void calculateScores(Scans& scans, TFMethod eMethod, float augmentFactor) {
    auto memes = &scans.getMemes();
    auto sequencesCount = &scans.getSequencesCount();

    auto memesIter = memes->begin();
    auto memesEnd = memes->end();
    double memesCount = (double)memes->size();
    // cout << "memes count: " << memesCount << endl;
    while (memesIter != memesEnd) {
        auto samplesIter = memesIter->second->getSamples().begin();
        auto samplesEnd = memesIter->second->getSamples().end();
        while (samplesIter != samplesEnd) {
            auto maxSequenceCount = samplesIter->second->getMaxSequenceCount();
            auto sequencesIter = samplesIter->second->getSequences().begin();
            auto sequencesEnd = samplesIter->second->getSequences().end();
            double score = 0;
            while (sequencesIter != sequencesEnd) {
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
                sequencesIter++;
            }
            samplesIter->second->setScore(score);
            samplesIter++;
        }
        memesIter++;
    }
}

void writeResults(Scans& scans, string& outputPath, string& bc) {
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
    
    auto memesIter = scans.getMemes().begin();
    auto memesEnd = scans.getMemes().end();
    int memeIndex = 0;
    while (memesIter != memesEnd) {
        hitsFile << "," << memesIter->first;
        scoreFile << "," << memesIter->first;

        auto samplesIter = memesIter->second->getSamples().begin();
        auto samplesEnd = memesIter->second->getSamples().end();
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
        
        auto memesIter = scans.getMemes().begin();
        auto memesEnd = scans.getMemes().end();
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

int main(int argc, char *argv[])
{
    cxxopts::Options options("td-idf", "Calculate td-idf based on hits");
    options.add_options()
        ("bc", "Biological condition", cxxopts::value<string>())
        ("sam2bc", "Samples to biological condition mapping", cxxopts::value<string>())
        ("scan", "Path to hits results directory", cxxopts::value<string>())
        ("output", "Path to results (CSVs) directory", cxxopts::value<string>())
        ("method", "TF method (boolean, terms, log, augmented), default is boolean", cxxopts::value<string>()->default_value("boolean"))
        ("factor", "Augment TF method factor (0-1), default is 0.5", cxxopts::value<float>()->default_value("0.5"))
        ("multiclass", "If labeling is multi-class or bc/other, default is false", cxxopts::value<bool>()->default_value("false"));

    auto result = options.parse(argc, argv);

    string bc = result["bc"].as<string>();
    string samples2bcPath = result["sam2bc"].as<string>();
    string scanPath = result["scan"].as<string>();
    string outputPath = result["output"].as<string>();
    string method = result["method"].as<string>();
    bool isMultiClass = result["multiclass"].as<bool>();
    float augmentFactor = result["factor"].as<float>();

    auto begin = chrono::steady_clock::now();

    auto eMethod = parseMethod(method);
    if (eMethod == TFMethod_NONE) {
        cout << "Invalid TF method!" << endl;
        return -1;
    }

    if (scanPath[scanPath.size() - 1] != '/') {
        scanPath += '/';
    }

    if (outputPath[outputPath.size() - 1] != '/') {
        outputPath += '/';
    }
    outputPath += bc + "/";
    mkdir(outputPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    Scans scans = loadScans(scanPath, samples2bcPath, bc, isMultiClass);
    calculateScores(scans, eMethod, augmentFactor);
    writeResults(scans, outputPath, bc);

    auto end = chrono::steady_clock::now();
    cout << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
}