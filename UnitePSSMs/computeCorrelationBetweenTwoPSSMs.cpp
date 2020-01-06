#include "computeCorrelationBetweenTwoPSSMs.h"
#include <vector>
#include <iostream>

// for pearson correlation
#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

void pushReleventPositionsToVectors(vector<double> & v1, vector<double> & v2, const PSSM & pssm1, const PSSM & pssm2, size_t p1, size_t p2) {
	for (size_t i = 0; i < pssm1.alphabetSize(); ++i) {
		v1.push_back(pssm1.getValue(p1, i));
		v2.push_back(pssm2.getValue(p2, i));
	}
}

double computeCorrelationBetweenTwoPSSMs(const MSA& inputMSA, const PSSM & pssm1, const PSSM & pssm2) {
	vector<double> v1, v2; // these two vectors will hold the relevant info in the matching PSSMs.
	size_t p1 = 0;
	size_t p2 = 0; // p1 and p2 are the current positions of the first and second sequence.
	double correlationRes = 0.0;
	for (size_t i = 0; i < inputMSA.getMSAlength(); ++i) {
		if (inputMSA.isMatch(i)) {
			pushReleventPositionsToVectors(v1, v2, pssm1, pssm2, p1, p2);
			p1++;
			p2++;
		}
		else if (inputMSA.isGap(0, i)) {// the first sequence has a gap
			p2++;
		}
		else if (inputMSA.isGap(1, i)) {
			p1++;
		}
		else {
			cerr << "you should not be here, function computeCorrelationBetweenTwoPSSMs" << endl;
		}
	}
	correlationRes = pearsoncoeff(v1, v2);
	return correlationRes;
}


// taken from: https://codepad.co/snippet/MbadrcBL
double sum(vector<double> a)
{
	double s = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		s += a[i];
	}
	return s;
}

double mean(vector<double> a)
{
	return sum(a) / a.size();
}

double sqsum(vector<double> a)
{
	double s = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		s += pow(a[i], 2);
	}
	return s;
}

double stdev(vector<double> nums)
{
	double N = nums.size();
	return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

vector<double> operator-(vector<double> a, double b)
{
	vector<double> retvect;
	for (size_t i = 0; i < a.size(); i++)
	{
		retvect.push_back(a[i] - b);
	}
	return retvect;
}

vector<double> operator*(vector<double> a, vector<double> b)
{
	vector<double> retvect;
	for (size_t i = 0; i < a.size(); i++)
	{
		retvect.push_back(a[i] * b[i]);
	}
	return retvect;
}

double pearsoncoeff(vector<double> X, vector<double> Y)
{
	return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}