#include "needleman_wunsch.h"


// this function is only used to test the implementation of the NW algorithm
// That is, it is used for debugging and QA when upgrading the source code.
void testNeedleman_wunsch() {
//int main() { 
	string testString1, testString2;
	vector<string> vstring1;
	
	cout << "Testing NW when both sequences are empty" << endl;
	testString1 = ""; 
	testString2 = "";
	needleman_wunsch nw1(2); // 2 for proteins.
	MSA msa1 = nw1.computeAffineNWForPair(testString1, testString2);
	msa1.printAligned();

	cout << "Testing NW when one sequence is empty" << endl;
	testString1 = "HAIM";
	testString2 = "";
	MSA msa2 = nw1.computeAffineNWForPair(testString1, testString2);
	msa2.printAligned();

	cout << "Testing NW when the second sequence is empty" << endl;
	testString1 = "";
	testString2 = "TAL";
	MSA msa3 = nw1.computeAffineNWForPair(testString1, testString2);
	msa3.printAligned();

	cout << "Testing NW for real sequences" << endl;
	//testString1 = "MLLALVCLLSCLLPSSEAKLYGRCELARVLHDFGLDGYRGYSLADWVCLAYFTSGFNAAALDYEADGSTNNGIFQINSRRWCSNLTPNVPNVCRMYCSDLLNPNLKDTVICAMKITQEPQGLGYWEAWRHHCQGKDLTEWVDGCDF"; // human lysozyme
	//testString2 = "MLNNRIFLLLSVFSMFSMCHIGAQIEVENKPVTEACLECLCEAMSGCNATKICVNGACGIFRITWTFWQDGGSLIGPGDGDAFTNCVNDPHCAADTIQNYMYKNGEDCNGDGKINCKDYGSIHKLGNLKCRDELPATFGTLFFKCLARKEQEEAEAQKKNST"; // drosophila lysozyme

	testString1 = "MLLALVCLLSCLLPSSEAKLYG"; // human lysozyme
	testString2 = "MLNNRIFLLLSVFSMFSMCH"; // drosophila lysozyme

	//testString1 = "MLLA"; // human lysozyme
	//testString2 = "MLN"; // drosophila lysozyme
	
	nw1.setGapExtensionScore(3);
	nw1.setGapOpenScore(7);
	MSA msa4 = nw1.computeAffineNWForPair(testString1, testString2);
	msa4.printAligned();
}


needleman_wunsch::needleman_wunsch(const vector<string> & seqArray,
	int similarity_mode,
	int match_score,
	int mismatch_score,
	int gap_open,
	int gap_extend) :
	_seqs(seqArray), _similarity_mode(similarity_mode), _match_score(match_score), _mismatch_score(mismatch_score), _gap_open(gap_open), _gap_extend(gap_extend)
{
	_simMat = sim_matrices(_similarity_mode, _match_score, _mismatch_score);
	_numberOfSequences = _seqs.size();
} 


needleman_wunsch::needleman_wunsch( //empty constructor, without sequences.
	int similarity_mode, // when this is 1, we align nucleotide. When 2, AA using Blosum62.
	int match_score,
	int mismatch_score,
	int gap_open,
	int gap_extend) :
	 _similarity_mode(similarity_mode), _match_score(match_score), _mismatch_score(mismatch_score), _gap_open(gap_open), _gap_extend(gap_extend)
{
	_simMat = sim_matrices(_similarity_mode, _match_score, _mismatch_score);
	_numberOfSequences = 0;
}

/*
needleman_wunsch::needleman_wunsch(const needleman_wunsch& other )//copy constructor
{
	_seqs = other._seqs;
	_numberOfSequences = other._numberOfSequences;
	_match_score = other._match_score;
	_mismatch_score = other._mismatch_score;
	_gap_open = other._gap_open;
	_gap_extend = other._gap_extend;
	_simMat = other._simMat;
	_similarity_mode = other._similarity_mode;
}*/

void needleman_wunsch::computeAllPairsMSAs(vector<MSA> & pairwiseMsasToFill)
{
	for(size_t i = 0; i < _numberOfSequences; i++)
	{
		for(size_t j = i + 1; j < _numberOfSequences; j++)
		{
			MSA curr_pairwise = computeAffineNWForPair(_seqs[i],_seqs[j]); //linear can be achieved if gap_open=gap_extend
			pairwiseMsasToFill.push_back(curr_pairwise);
		}
	}
}


MSA needleman_wunsch::computeAffineNWForPair(const string & A, const string & B) {
	size_t m = A.length();
	size_t n = B.length();

	// deal with edge cases where one or more of the sequences is of size zero:
	if((A.length() == 0) && (B.length() == 0)) {
		vector<string> aligned_seqs;
		aligned_seqs.push_back("");
		aligned_seqs.push_back("");
		MSA curr_pairwise = MSA(aligned_seqs, INT_MIN);
		return curr_pairwise;
	}

	if(A.length() == 0) {
		string retA = "";
		retA.insert(0, B.length(), '-'); // if A is empty we enter B.length() instances of '-' to the aligned sequence.
		vector<string> aligned_seqs;
		aligned_seqs.push_back(retA);
		aligned_seqs.push_back(B);
		int aln_score = _gap_open + (B.length() - 1)*_gap_extend;
		MSA curr_pairwise = MSA(aligned_seqs, aln_score);
		return curr_pairwise;
	}

	if(B.length() == 0) {
		string retB = "";
		retB.insert(0, A.length(), '-'); 
		vector<string> aligned_seqs;
		aligned_seqs.push_back(A);
		aligned_seqs.push_back(retB);
		int aln_score = _gap_open + (A.length() - 1)*_gap_extend;
		MSA curr_pairwise = MSA(aligned_seqs, aln_score);
		return curr_pairwise;
	}
	// end deal with edge cases where one or more of the sequences is of size zero
	// This computation is based on the book by Durbin (Biological sequence analysis), pg 29 
	
	// The M matrix. M[i][j] is the best score when the first sequence is aligned till position i,
	// the second sequence is aligned till position j, and position i is aligned to position j.
	vector<vector<int>> M_mat(m + 1, vector<int>(n + 1)); 
	
	// The I matrix. I[i][j] is the best score when the first sequence is aligned till position i,
	// the second sequence is aligned till position j, and position i in the first sequence is
	// aligned to a gap.
	vector<vector<int>> I_mat(m + 1, vector<int>(n + 1));
	
	vector<vector<int>> H_mat(m+1, vector<int>(n+1));
	
	
	vector<vector<int>> J_mat(m+1, vector<int>(n+1));

	//matrices initialization

	// Initialization of the M matrix. In M, we assume that position i is aligned to position j.
	// It cannot be that position 2 of the first sequence is aligned with a match to position 0
	// of the second sequence. So, we set M[0,i] and M[j,0] to be minus infinity.

	// The I matrix captures the cases in which the upper sequence is extended
	// For example    AAKC
	// Aligned with   A_T_
	// I[0][0] is impossible and hence it is set to minus infinity.
	// I[0][5] is also impossible and hence it is set to minus inf as well.
	// I[5][0] suggests that there are five gaps in the second sequence, so it
	// is set to gap openning + the cost of four gaps.

	// The matrix H stores the max of the matrices M,I,J.
	for(size_t j=1; j<= B.length(); j++) {
		M_mat[0][j] = -(INT_MAX / 2);
		I_mat[0][j] = -(INT_MAX / 2);
		J_mat[0][j] = -_gap_open - _gap_extend *(j - 1);
		H_mat[0][j] = - _gap_open - _gap_extend *(j - 1);	
	}

	for(size_t i=1; i<=m; i++) {
		M_mat[i][0] = -(INT_MAX / 2);
		I_mat[i][0] = -_gap_open - _gap_extend *(i - 1);
		J_mat[i][0] = -(INT_MAX / 2);
		H_mat[i][0] = - _gap_open - _gap_extend *(i - 1);	
	}

	H_mat[0][0] = 0;
	M_mat[0][0] = 0; // this is a special case of M, which is ok (zero can be aligned to zero).
	I_mat[0][0] = - (INT_MAX /2);
	J_mat[0][0] = - (INT_MAX /2);
	// end initialization

	//score computation - dynamic programming
	for(size_t i=1; i<=m; i++) {
		for(size_t j=1; j<=n; j++) {
			int S = _simMat.get_pair_score(A[i-1],B[j-1]);
			int M_diag = M_mat[i-1][j-1] + S;
			int I_diag = I_mat[i-1][j-1] + S;
			int J_diag = J_mat[i-1][j-1] + S;
			M_mat[i][j] = max(M_diag, max(I_diag, J_diag));
			
			// here you insert a gap in the top sequence
			int M_up = M_mat[i-1][j];
			int I_up = I_mat[i-1][j];
			I_mat[i][j] = max((M_up - _gap_open),(I_up - _gap_extend));

			int M_left = M_mat[i][j-1];
			int J_left = J_mat[i][j-1];
			J_mat[i][j] = max((M_left - _gap_open),(J_left - _gap_extend));

			H_mat[i][j] = max(M_mat[i][j], max(I_mat[i][j], J_mat[i][j]));
		}
	}
	// end score computation - dynamic programming

	// traceback
	string retA, retB;
    stack<char> SA, SB;

    int ii = m;
	int jj = n;

    while (ii != 0 || jj != 0)
    {
        if (ii == 0)
        {
            SA.push('-');
            SB.push(B[jj-1]);
            jj--;
        }
        else if (jj == 0)
        {
            SA.push(A[ii-1]);
            SB.push('-');
            ii--;
        }
        else
        {
			int curr_H_val = H_mat[ii][jj];
			if(curr_H_val == M_mat[ii][jj]) {
				//go diag
				SA.push(A[ii-1]); // A [ii-1] is the ii'th position in sequence A.
                SB.push(B[jj-1]);
                ii--;
				jj--;
			}
			else if(curr_H_val == I_mat[ii][jj])
			{
				//take from A
				SA.push(A[ii-1]);
                SB.push('-');
                ii--;
			}
			else
			{
				//take from B
				SA.push('-');
                SB.push(B[jj-1]);
                jj--;
			}
        }
    }

	while (!SA.empty())
    {
        retA += SA.top();
        retB += SB.top();
        SA.pop();
        SB.pop();
    }
	vector<string> aligned_seqs;
	aligned_seqs.push_back(retA);
	aligned_seqs.push_back(retB);
	MSA curr_pairwise(aligned_seqs, H_mat[m][n]);
	//cout << "The score of the best alignment = " << H_mat[m][n] << endl;
    return curr_pairwise;

}