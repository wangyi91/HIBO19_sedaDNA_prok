#!/usr/bin/env julia

using Arrow
using DataFrames
using GZip
using CSVFiles
using Pipe: @pipe

include("env.jl") # WORK_DIR
include("./InitialExploration/_initial_load.jl")
include("./MicrobeProfiling/_get_annotated_df.jl") # to use the tax_at_rank function

pip="amaw"
alg = ".lca";
add="_ANI92";
tag = "$pip$alg$add"

# ------------------ ONLY NEED TO RUN THIS ONCE ----------------------- #
write_arrow(tag)
#-----------------------------------------------------------------------#


# join metadata; filter based on N_reads and tax_rank
nreads=50; ftax=""; rank="species";
ftax = ["Bacteria","Archaea"]

## data of samples
(x,y) = load_arrow_as_df(tag, ftax, rank, nreads)

# read acc2taxid file
using CSVFiles
acc2taxid = DataFrame(load(File(format"TSV", WORK_DIR*"vanilla-organelles-virus_7nov2022/pkg/taxonomy/acc2taxid.map.gz")))

# load mapping results
tb = CSV.read(WORK_DIR*"HIBO_shotgun/EnvironmenTracker/archive_results/results_latest_amaw$add/taxonomic-profiling/tp-mapping-filtered.summary.tsv.gz", DataFrame, delim='\t');
tb1 = @pipe tb |> leftjoin(_, acc2taxid[!,[:taxid, :accession]], on=:reference=>:accession)

# for microbe analysis: remove plastid and mito data to view only prokaryotes
filter!(:reference => r-> !contains(r, "_mito"),tb);
filter!(:reference => r-> !contains(r, "_plas"),tb);


# create a data file containing both samples and controls
z = vcat(x, y)
z.tax_name = [n[4:end] for n in z.tax_name]
leftjoin!(z,tb1, on=[:tax_id=>:taxid, :Label=>:label])

x.tax_name = [n[4:end] for n in x.tax_name]
leftjoin!(x,tb1, on=[:tax_id=>:taxid, :Label=>:label])


using JLD2
save_object("./InitialExploration/data/$tag.$rank.jld2",x)
save_object("./InitialExploration/data/$tag.$rank.samplecontrol.jld2",z)

