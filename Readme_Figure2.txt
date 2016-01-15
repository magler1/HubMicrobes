# To run the scripts, one needs to have installed:
# 1. R
# 2. Perl
# 3. The perl package: Statistics-R-0.33 and required dependencies (including IPC-Run-0.94 and Regexp-Common-2013031301)
# Also, the directions are set up for running on a unix-based environment, so need to be adapted otherwise

# First source the required perl packages if you need to do so, for example:
# export PERL5LIB=$PERL5LIB:/Users/agler/packages/Statistics-R-0.33/lib/:/Users/agler/packages/IPC-Run-0.94/lib/:/Users/agler/packages/Regexp-Common-2013031301/lib/

# Run the script to calculate the raw correlations
perl Find_OTU_Correlations_NoFilt.pl otu_correlations_list_all.txt

# Concatenate the results from all of the comparisons together into master node and edge files

mkdir ./Correlation_Resubmission_Test_NoFilt/All_NodeTables/
mkdir ./Correlation_Resubmission_Test_NoFilt/All_EdgeTables/

ls ./Correlation_Resubmission_Test_NoFilt/*vs*/Node*4.txt | xargs cat > ./Correlation_Resubmission_Test_NoFilt/All_NodeTables/Master_NodeTable4.txt
ls ./Correlation_Resubmission_Test_NoFilt/*vs*/Node*5.txt | xargs cat > ./Correlation_Resubmission_Test_NoFilt/All_NodeTables/Master_NodeTable5.txt
ls ./Correlation_Resubmission_Test_NoFilt/*vs*/Node*6.txt | xargs cat > ./Correlation_Resubmission_Test_NoFilt/All_NodeTables/Master_NodeTable6.txt
ls ./Correlation_Resubmission_Test_NoFilt/*vs*/Edge*4.txt | xargs cat > ./Correlation_Resubmission_Test_NoFilt/All_EdgeTables/Master_EdgeTable4.txt
ls ./Correlation_Resubmission_Test_NoFilt/*vs*/Edge*5.txt | xargs cat > ./Correlation_Resubmission_Test_NoFilt/All_EdgeTables/Master_EdgeTable5.txt
ls ./Correlation_Resubmission_Test_NoFilt/*vs*/Edge*6.txt | xargs cat > ./Correlation_Resubmission_Test_NoFilt/All_EdgeTables/Master_EdgeTable6.txt

# Cleanup the node and edge tables to filter as necessary and to get rid of self-correlations (use this as many times as necessary, modifying the output
# file name to generate filtered networks at various cutoff strengths. The settings "as-is" will generate a network table similar to Figure 2

perl CleanUp_Tables_Simple_NoSelfEdges.pl

# Note that the previous script might have left the column names in the edge and node tables in the wrong row. One should move this row to the top.
# Load the edge and node tables into cytoscape and calculate statistics on the network. These can then be output as a file with properties of each node.
# Cytoscape settings for figure2 (use the L6, genus-level, network):

# Edge Property|Parameter|Cont/Discrete|Value-Setting|Value-Setting
# Type| Slope|Discrete|Neg-Dashed|Pos-Solid
# Color|Slope|Discrete|Neg-Black|Pos-Orange
# Transparency|stm_rsqadj|Continuous|3-30|7.0274-15
# Width|stm_rsqadj|Continuous|3-0|7.0274-15

# Node Property|Parameter|Cont/Discrete|Value-Setting|Value-Setting|Value-Setting
# Size|Degree|Continuous|0-10|53-50|
# Color|Kingdom|Discrete|Bacteria-Blue|Chromalveolata-Green|Fungi-Red
# Border|EndoEpi|Discrete|Endo-Black|Epi-None
