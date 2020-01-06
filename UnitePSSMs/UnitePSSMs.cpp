#include "PSSM.h"
#include "fromFileToPSSM_array.h"
#include "needleman_wunsch.h"
#include "MSA.h"
#include "computeCorrelationBetweenTwoPSSMs.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cfloat>
using namespace std;

int mainPipeline(int argc, char *argv[]);

int main(int argc, char *argv[]) {
	mainPipeline(argc,argv);
}

using namespace std;

vector<vector<PSSM> > cluster_Motifs(vector<PSSM> pssmVec, double correlationCutof, int Aln_Score_Cutoff, string log_file);
void print_CSV_list_of_motif_to_unite(string FileName, vector<vector<PSSM>> & PSSM_to_unite);
void sortPSSMarrayAccordingToNsite(vector<PSSM>& PSSM_array);
//struct less_than_PSSM_Nsite
//{
//	inline bool operator() (const PSSM& pssm1, const PSSM& pssm2)
//	{
//		return (pssm1.get_Nsite() < pssm2.get_Nsite());
//	}
//};


void sortPSSMarrayAccordingToNsite(vector<PSSM>& PSSM_array) {
	//	sort(PSSM_array.begin(), PSSM_array.end(), less_than_PSSM_Nsite());
	sort(PSSM_array.begin(), PSSM_array.end(), [](const PSSM& lhs, const PSSM& rhs)
	{
		return lhs.get_Nsite() > rhs.get_Nsite();
	});
}

int mainPipeline(int argc, char *argv[])
{
	string In_PSSM_FileName = "";
	string Out_FileName = "";
	int Aln_Score_Cutoff = INT_MIN;
	double PCC_Cutoff = DBL_MAX;
	
	int num_required_params = 4;

	if (argc < (num_required_params * 2)) // each with its flag
	{
		// Check the value of argc. If not enough parameters have been passed, inform user and exit.
		std::cout << "Usage is "<<argv[0]<<" -pssm <PSSMs_in_MAST_Format> -out <out_file_for_list_of_PSSMs_to_unite> -aln_cutoff <cutoff for pairwise alignment score to unite> -pcc_cutoff <minimal PCC R to unite>\n"; // Inform the user of how to use the program
		std::cin.get();
		exit(1);
	}
	else
	{
		// if we got enough parameters...
		// char* myFile, myPath, myOutPath;
		cout << argv[0];
		for (int i = 1; i < (num_required_params * 2); i += 2) // iterate over flags
		{
			/* We will iterate over argv[] to get the parameters stored inside.
			* Note that we're starting on 1 because we don't need to know the
			* path of the program, which is stored in argv[0] */
			if (string(argv[i]) == "-pssm")
			{
				// We know the next argument *should* be the filename:
				In_PSSM_FileName = string(argv[i + 1]);
				cout << " -pssm " << In_PSSM_FileName;
			}
			else if (string(argv[i]) == "-out")
			{
				Out_FileName = string(argv[i + 1]);
				cout << " -out " << Out_FileName ;
			}
			else if (string(argv[i]) == "-aln_cutoff")
			{
				Aln_Score_Cutoff = atoi(argv[i + 1]);
				cout << " -aln_cutoff " << Aln_Score_Cutoff;
			}
			else if (string(argv[i]) == "-pcc_cutoff")
			{
				PCC_Cutoff = atof(argv[i + 1]);
				cout << " -pcc_cutoff " << PCC_Cutoff;
			}
			else
			{
				std::cout << "Not enough or invalid arguments, please try again.\n";
				exit(1);
			}
      cout<<endl;
		}
	}
	string log_file = Out_FileName + ".log.txt";
	vector<PSSM> PSSM_array;
//	cout << "Hello world!" << endl;
	map<string, int>  PSSM_Name_To_ArrayIndex; // map between the location of the pssm in the array and its name.
//	In_PSSM_FileName = "C:/Data/Dropbox/haim/motifs/UnitePSSMs/Her_motifs.txt";
	readFileToPSSM_array(In_PSSM_FileName, PSSM_array, PSSM_Name_To_ArrayIndex);
	sortPSSMarrayAccordingToNsite(PSSM_array);
//	vector<string> seq_to_align;
//	string seqA = "A";
//	string seqB = "C";
//	seq_to_align.push_back(seqA);
//	seq_to_align.push_back(seqB);
//	needleman_wunsch NW_obj_test(seq_to_align, 2);
//	MSA test = NW_obj_test.computeAffineNWForPair(seqA, seqB);

	// testing random peptide values
	// PSSM testPSSM = PSSM_array[0];
	// testPSSM.giveScorePercentile(3, 0.95, 100);
	// return 0;

	// the code to unite PSSMs:
	vector<vector<PSSM>> PSSM_to_unite = cluster_Motifs(PSSM_array, PCC_Cutoff,Aln_Score_Cutoff,log_file);
	print_CSV_list_of_motif_to_unite(Out_FileName, PSSM_to_unite);
//	cout << "Hello world" << endl;
	return 1;
}
vector<vector<PSSM> > cluster_Motifs(vector<PSSM> pssmVec, double correlationCutof, int Aln_Score_Cutoff, string log_file) {// note that we copy the input pssm 
	// we returned a vector of PSSM vectors. Each entry is a list of PSSMs we should unite.
	//This function does not actually unit the PSSMs.
	vector<vector<PSSM> > res;
	ofstream LOG_FILE;
	LOG_FILE.open(log_file);
	while (pssmVec.size()>0) {
		LOG_FILE << "== " << pssmVec[0].get_PSSM_name() << " is a new cluster head" << endl;
		vector<PSSM> tmp; //this will hold each cluster as we build it
		for (size_t i = 1; i < pssmVec.size(); ++i) {
			vector<string> seq_to_align;
			seq_to_align.push_back(pssmVec[0].getConsensus());
			seq_to_align.push_back(pssmVec[i].getConsensus());
			needleman_wunsch NW_obj_test(seq_to_align, 2);
			MSA alignedMSA = NW_obj_test.computeAffineNWForPair(pssmVec[0].getConsensus(), pssmVec[i].getConsensus());
			if (alignedMSA.getAlignmentScore() > Aln_Score_Cutoff)
			{
				double PSSM_Correlation = computeCorrelationBetweenTwoPSSMs(alignedMSA, pssmVec[0], pssmVec[i]);
				if (PSSM_Correlation > correlationCutof) {
					LOG_FILE << "Added PSSM: " << pssmVec[i].get_PSSM_name() << " with Pearson correlation of " << PSSM_Correlation << " with alignment score " << alignedMSA.getAlignmentScore() << " and alignment " << endl << alignedMSA.getAlignedSeqs()[0] << endl << alignedMSA.getAlignedSeqs()[1] << endl;
					tmp.push_back(pssmVec[i]);
					pssmVec.erase(pssmVec.begin() + i); // remove this PSSM from PSSM list left to unite
					i--; // because the element was removed, we want i to go back...
				}
			}
		}
		tmp.push_back(pssmVec[0]);
		pssmVec.erase(pssmVec.begin());
		res.push_back(tmp);
		LOG_FILE << "Total motives in cluster: " << tmp.size()<< endl;
		LOG_FILE << "===================================================================================" << endl<<endl;
		cout << "== PSSMs left to cluster: " << pssmVec.size()<<endl;
	}
	return res;
}
void print_CSV_list_of_motif_to_unite(string FileName, vector<vector<PSSM>> & PSSM_to_unite)
{
	ofstream CSV_file;
	CSV_file.open(FileName);
	for (size_t i = 0; i < PSSM_to_unite.size(); ++i)
	{
		for (size_t j = 0; j < PSSM_to_unite[i].size(); j++)
		{
			CSV_file << PSSM_to_unite[i][j].get_PSSM_name();
			if (j != PSSM_to_unite[i].size()-1) { CSV_file << "," ; }
			else { CSV_file << endl; }
		}
	}
	CSV_file.close();
}