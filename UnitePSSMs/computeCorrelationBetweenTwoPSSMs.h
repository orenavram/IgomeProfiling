#ifndef _COMPUTE_CORR_TWO_PSSM
#define _COMPUTE_CORR_TWO_PSSM
#include "MSA.h"
#include "PSSM.h"
#include <string>
using namespace std; 

double sum(vector<double> a);
double mean(vector<double> a);
double sqsum(vector<double> a);
double stdev(vector<double> nums);
vector<double> operator-(vector<double> a, double b);
vector<double> operator*(vector<double> a, vector<double> b);
double pearsoncoeff(vector<double> X, vector<double> Y);

double computeCorrelationBetweenTwoPSSMs(const MSA& inputMSA, const PSSM & pssm1, const PSSM & pssm2);




#endif

