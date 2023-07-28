using DrWatson
@quickactivate "FindrTutorials"

using DataFrames
using CSV
using Arrow
using Printf



#===================
Read data

Note: because julia is column-major, and most of our analyses work on gene expression profiles across samples, expression data is stored with genes as columns and samples as rows
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
Arrow.write(oname,expr)

oname = datadir("exp_pro","findr-data-geuvadis", "dm.arrow")
Arrow.write(oname,mirna)

oname = datadir("exp_pro","findr-data-geuvadis", "dgt.arrow")
Arrow.write(oname,geno)

oname = datadir("exp_pro","findr-data-geuvadis", "dgm.arrow")
Arrow.write(oname,geno_mirna)
