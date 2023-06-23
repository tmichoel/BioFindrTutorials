#############################################
#           Download instructions           #
#############################################

# 1. Create a directory data/exp_raw/ecoli-grn

# 2. Download expression data from
# - URL: http://m3d.mssm.edu/norm/
# - Download: E_coli_v4_Build_6.tar.gz
# - Unzip and extract (tar -xvf) the tar file, which will create a sub-dir E_coli_v4_Build_6
# - We will use the file avg_E_coli_v4_Build_6_exps466probes4297.tab


# 3. Download ground-truth GRN from
# - URL: https://regulondb.ccg.unam.mx/menu/download/datasets/index.jsp
# - Download: Regulatory Network Interactions -> TF-gene interactions, or diretly from this URL
# https://regulondb.ccg.unam.mx/menu/download/datasets/files/NetWorkTFGene.txt


##################################################################
#           Activate the environment and load packages           #
##################################################################

using DrWatson
@quickactivate "FindrTutorials"

using DataFrames
using CSV
using Arrow

#############################################
#           Gene expression data            #
#############################################

# Data file 
fexpr = datadir("exp_raw","ecoli-grn","E_coli_v4_Build_6", "avg_E_coli_v4_Build_6_exps466probes4297.tab");

# Read into dataframe, truncate gene names at first "_", and remove the experiment IDs (1st column)
dfexpr = permutedims(DataFrame(CSV.File(fexpr)),1);
rename!(x -> split(x,"_")[1], dfexpr);
select!(dfexpr, Not(1));

###########################################
#           RegulonDB GRN data            #
###########################################

# Data file 
fgrn = datadir("exp_raw","ecoli-grn","NetWorkTFGene.txt");

# Read into dataframe, ignore commented lines
dfgrn = DataFrame(CSV.File(fgrn; comment="#"))

# Keep only columns with regulator and regulated gene name, and confidencelevel
select!(dfgrn, :"3)RegulatorGeneName", :"5)regulatedName", :"7)confidenceLevel")

# Rename columns
rename!(dfgrn, [:TF, :Target, :Confidence])

#######################################################
#           Align expression and GRN data            #
#######################################################

# Keep only GRN rows where TF and target both have expression data
filter!(x -> ∈(x.TF, names(dfexpr)) & ∈(x.Target, names(dfexpr)), dfgrn)

# Keep only genes in expression data that are in GRN
ingrn = map(x -> ∈(x,  union(dfgrn.TF,dfgrn.Target)) , names(dfexpr))
select(dfexpr, ingrn)


#####################################
#       Save processed data         #
#####################################

fexpr_pro = datadir("exp_pro","ecoli-grn","avg_E_coli_v4_Build_6_exps466probes4297_filtered.arrow");
Arrow.write(fexpr_pro,dfexpr)

fgrn_pro = datadir("exp_pro","ecoli-grn","RegulonDB_v11.1_NetworkTFGene.arrow");
Arrow.write(fgrn_pro,dfgrn)