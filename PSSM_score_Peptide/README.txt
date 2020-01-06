# To compile
module load gcc/gcc620
/groups/pupko/orenavr2/igomeProfilingPipeline/src/PSSM_score_Peptide
g++ *.cpp -std=c++11 -O3 -o PSSM_score_Peptide

# To run
1. To calculate the cutoffs for each PSSM
./PSSM_score_Peptide -pssm example/example_Motifs_File.meme_format -pssm_cutoffs example/example_PSSM_Scores_Cutoffs.txt -CalcPSSM_Cutoff

2. Calculate PSSM Pvalues [using 100 randomized PSSMs]
./PSSM_score_Peptide -pssm example/example_Motifs_File.meme_format -pssm_cutoffs example/example_PSSM_Scores_Cutoffs.txt -seq example/example_SeqFile.fas -out example/example_PSSMs_Pvals.txt -NrandPSSM 100 -CalcPSSM_Pval 

3. find hits using cutoffs
./PSSM_score_Peptide -pssm example/example_Motifs_File.meme_format -pssm_cutoffs example/example_PSSM_Scores_Cutoffs.txt -seq example/example_SeqFile.fas -out example/example_PSSM_hits.txt -CalcPSSM_Hits 
