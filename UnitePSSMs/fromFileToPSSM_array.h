#ifndef FROM_FILE_TO_PSSM_____H
#define FROM_FILE_TO_PSSM_____H

#include "PSSM.h"
#include <string>
#include <vector>
#include <map>
using namespace std;

void readFileToPSSM_array(const string fileName, vector<PSSM> & PSSM_array, map<string, int> & PSSM_Name_To_ArrayIndex);






#endif