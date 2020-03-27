#include <string>
#include <ostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>

#include "cxxopts.hpp"
#include "types.hpp"
#include "memeSample.hpp"
#include "meme.hpp"
#include "scans.hpp"
#include "utils.hpp"

Scans loadScans(string scanPath, string samples2bcPath, string& bc, bool isMultiClass);
void calculateScores(Scans& scans, TFMethod eMethod, float augmentFactor);
void writeResults(Scans& scans, string& outputPath, string& bc);

int main(int argc, char *argv[])
{
    cxxopts::Options options("td-idf", "Calculate td-idf based on hits");
    options.add_options()
        ("bc", "Biological condition", cxxopts::value<string>())
        ("sam2bc", "Samples to biological condition mapping", cxxopts::value<string>())
        ("scan", "Path to hits results directory", cxxopts::value<string>())
        ("output", "Path to results (CSVs) directory", cxxopts::value<string>())
        ("done", "Path to done file", cxxopts::value<string>())
        ("method", "TF method (boolean, terms, log, augmented), default is boolean", cxxopts::value<string>()->default_value("boolean"))
        ("factor", "Augment TF method factor (0-1), default is 0.5", cxxopts::value<float>()->default_value("0.5"))
        ("multiclass", "If labeling is multi-class or bc/other, default is false", cxxopts::value<bool>()->default_value("false"))
        ("v,verbose", "Verbose output", cxxopts::value<string>()->default_value("false"));

    auto result = options.parse(argc, argv);

    string bc = result["bc"].as<string>();
    string samples2bcPath = result["sam2bc"].as<string>();
    string scanPath = result["scan"].as<string>();
    string outputPath = result["output"].as<string>();
    string donePath = result["done"].as<string>();
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

    ofstream doneFile(donePath);
    doneFile << "Calculated TF-IDF for " << bc << " biological condition" << endl;
    doneFile << "Completed in " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
}
