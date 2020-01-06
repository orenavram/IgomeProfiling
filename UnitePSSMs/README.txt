# To compile
module load gcc/gcc620
cd /groups/pupko/orenavr2/igomeProfilingPipeline/src/UnitePSSMs
g++ *.cpp -std=c++11 -o UnitePSSMs

# To run
./UnitePSSMs -pssm <PSSMs_in_MAST_Format> -out <out_file_for_list_of_PSSMs_to_unite> -aln_cutoff <cutoff for pairwise alignment score to unite> -pcc_cutoff <minimal PCC R to unite>

