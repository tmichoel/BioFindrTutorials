using DrWatson
@quickactivate "BioFindrTutorials"

using DataFrames
using CSV
using Arrow
using Printf
using Statistics


#===================
Read data

Note: because julia is column-major, and most of our analyses work on gene expression profiles across samples, expression data is stored with genes as columns and samples as rows.
====================#

# Read mrna expression data from exp_raw folder
fname = datadir("exp_raw","findr-data-geuvadis","dt.csv")
expr = permutedims(DataFrame(CSV.File(fname)),1)

# Read mirna expression data from exp_raw folder
fname = datadir("exp_raw","findr-data-geuvadis","dm.csv")
mirna = permutedims(DataFrame(CSV.File(fname)),1)

# Read mrna genotypes from exp_raw folder
fname = datadir("exp_raw","findr-data-geuvadis","dgt.csv")
geno = permutedims(DataFrame(CSV.File(fname)),1;makeunique=true)

# Read mirna genotypes from exp_raw folder
fname = datadir("exp_raw","findr-data-geuvadis","dgm.csv")
geno_mirna = permutedims(DataFrame(CSV.File(fname)),1;makeunique=true)

# Verify that all rows are aligned
sum(expr.Column1 .== mirna.Column1) / nrow(expr)
sum(expr.Column1 .== geno.Column1) / nrow(expr)
sum(expr.Column1 .== geno_mirna.Column1) / nrow(expr)

# We can drop the sample names
select!(expr, Not(:Column1))
select!(mirna, Not(:Column1))
select!(geno, Not(:Column1))
select!(geno_mirna, Not(:Column1))

# Read conversion.txt file to map Ensembl IDs to gene names
fname = datadir("exp_raw","findr-data-geuvadis","conversion.txt")
gmap = DataFrame(CSV.File(fname))
rename!(gmap, [:"Ensembl Gene ID" => :EnsemblID, :"Associated Gene Name" => :GeneName])

enames = DataFrame(:EnsemblID => map(x -> split(x,".")[1], names(expr)))
leftjoin!(enames, gmap, on = :EnsemblID)

# Save data to exp_pro folder
oname = datadir("exp_pro","findr-data-geuvadis", "dt.arrow")
Arrow.write(oname,dt)

oname = datadir("exp_pro","findr-data-geuvadis", "dm.arrow")
Arrow.write(oname,mirna)

oname = datadir("exp_pro","findr-data-geuvadis", "dgt.arrow")
Arrow.write(oname,dgt)

oname = datadir("exp_pro","findr-data-geuvadis", "dgm.arrow")
Arrow.write(oname,geno_mirna)

#================
Rename gene names to HGNC symbols
=================#

# Load data
dt = DataFrame(Arrow.Table(datadir("exp_pro","findr-data-geuvadis", "dt.arrow")))
dgt = DataFrame(Arrow.Table(datadir("exp_pro","findr-data-geuvadis", "dgt.arrow")))

# Load conversion
dconv = DataFrame(CSV.File(datadir("exp_raw","findr-data-geuvadis","HGNC_ensembl_conversion.txt")))

# create a dictionary to map ensembl ids to approved symbols
map_names = Dict(zip(dconv[!,:"Ensembl gene ID"], dconv[!,:"Approved symbol"]))

# rename the columns
new_names = [get(map_names, x, x) for x in names(dt)]
rename!(dt, new_names)

# find the indices of the genes that were not renamed
idx = findall(x -> !∈(x,dconv.:"Approved symbol"), new_names)

# create the matching indices for dgt
idx_dgt = idx[idx .<= ncol(dgt)]

# remove the genes that were not renamed in dt and dgt
select!(dt, Not(idx))
select!(dgt, Not(idx_dgt))

# check the correlation between aligned columns
cc = [cor(dt[:,i], dgt[:,i]) for i in 1:ncol(dgt)]

# keep only the columns that have a correlation greater than 0.2
idx = findall(x -> abs(x) < 0.2, cc)
select!(dt, Not(idx))
select!(dgt, Not(idx))

#========================
Lingfei's conversion and gold standard data
=========================#

# Load conversion
dconv = DataFrame(CSV.File(datadir("exp_raw","findr-data-geuvadis","conversion.txt")))

# Rename gene names
map_names = Dict(zip(dconv[!,:"Ensembl Gene ID"], dconv[!,:"Associated Gene Name"]))
map_names_rev = Dict(zip(dconv[!,:"Associated Gene Name"], dconv[!,:"Ensembl Gene ID"]))

# split the ensembl id to remove the version number
new_names = map(x -> split(x, ".")[1], names(dt))
rename!(dt, new_names)

# replace ensembl id by gene name if possible
new_names = [get(map_names, x, x) for x in new_names]

# read gold standard data files
encode = DataFrame(CSV.File(datadir("exp_raw", "findr-data-geuvadis", "gold", "encode.csv"), header=false))
sirna = DataFrame(CSV.File(datadir("exp_raw", "findr-data-geuvadis", "gold", "sirna.csv")))
mirlab = DataFrame(CSV.File(datadir("exp_raw", "findr-data-geuvadis", "gold", "mirlab.csv")))

# get enembl ids from sirna names
TF_sirna = map(x -> split(x, "_")[2], names(sirna)[4:end])
intersect!(TF_sirna, names(dt))

# get encode TFs
TF_encode_names = unique(encode.Column1)
TF_encode = [get(map_names_rev, x, x) for x in TF_encode_names]
intersect!(TF_encode, names(dt))

# all TFs
TF = union(TF_sirna, TF_encode)

# save TFs to file
oname = datadir("exp_pro","findr-data-geuvadis", "TFs.csv")
CSV.write(oname, DataFrame(TF=TF))

# mircoRNAs in mirlab
miRNAs = intersect(unique(mirlab.miRNA), names(dm))

# save miRNAs to file
oname = datadir("exp_pro","findr-data-geuvadis", "miRNAs.csv")
CSV.write(oname, DataFrame(miRNA=miRNAs))