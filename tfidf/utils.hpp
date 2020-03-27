#pragma once
#include <string>
#include "types.hpp"

using namespace std;
void tolower(string& str);
TFMethod parseMethod(string method);
SamplesBC loadSamplesToBC(string samples2bcPath);