#!/usr/bin/env julia
using DataFrames, JLD2, Pipe, CSV
using FlashWeave
using Graphs#, SimpleWeightedGraphs
using Compose, GraphPlot, Colors
import Cairo, Fontconfig

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw";add="_ANI92";
alg=".lca";
tag="$pip$alg$add"
rank="species"


# Step 1: load data_training by running ./InitialExploration/vcat.jl 
otudata=load_object("./deContamination/data/$tag.$rank.jld2")

# Step 2: filter out singletons, convert data to otu table and save as csv
keep = @pipe otudata |> combine(groupby(_, :tax_path), nrow=>"N_occur") |> filter(:N_occur=>o->o>=2, _)

otu = @pipe otudata |> filter(:tax_path => p->p in keep.tax_path, _) |> filter(:N_reads => n->n>50, _) |>
        unstack(_, :Label, :tax_path, :N_reads, fill=0)

data_path = "./NetworkAnalysis/data/otu_$tag.csv"
CSV.write(data_path, otu[!,2:end])

    ##### Additionally, create 2 extra otu tables with sp in at least 3 and 5 samples ########
    keep3 = @pipe otudata |> combine(groupby(_, :tax_path), nrow=>"N_occur") |> filter(:N_occur=>o->o>=3, _)
    keep5 = @pipe otudata |> combine(groupby(_, :tax_path), nrow=>"N_occur") |> filter(:N_occur=>o->o>=5, _)

    otu3 = @pipe otudata |> filter(:tax_path => p->p in keep3.tax_path, _) |>
        DataFrames.transform(_, :N_reads=> ByRow(n->sqrt(n))=>:sqrtN_reads) |>
        unstack(_, :Label, :tax_path, :N_reads, fill=0) 
    otu5 = @pipe otudata |> filter(:tax_path => p->p in keep5.tax_path, _) |>
        DataFrames.transform(_, :N_reads=> ByRow(n->sqrt(n))=>:sqrtN_reads) |>
        unstack(_, :Label, :tax_path, :N_reads, fill=0) 

    data_path3="./NetworkAnalysis/data/otu3_$tag.csv";
    data_path5="./NetworkAnalysis/data/otu5_$tag.csv";
    CSV.write(data_path3, otu3[!,2:end])
    CSV.write(data_path5, otu5[!,2:end])
    ##########################################################################################


# Step 3: compile metadata and save as csv
    meta = compile_metadata()
    save_object("./NetworkAnalysis/data/meta.jld2", meta)
# easy access
meta = load_object("./NetworkAnalysis/data/meta.jld2")

## final filtering
md = @pipe meta |> DataFrames.transform(_, :yrBP=>ByRow(y->y*(-1))=>:yrBP) |> # for better interpretation of connection with time 
        DataFrames.transform!(_, :Middle_depth => ByRow(d -> split(string(d),".")[1]*"cm") => :d) |> # add depth as a categorical MV
        select(_,Not([:Middle_depth,:sigma])); md.S[8]=142; md.bSi[33]=md.bSi[34]=1.5;# manually replace missing values

## normalise MVs
md[!,[:S,:bSi,:TOC,:TIC,:yrBP]] = Float64.(md[:,[:S,:bSi,:TOC,:TIC,:yrBP]]); # change column type
using StatsBase
md[!,[:S,:bSi,:TOC,:TIC,:yrBP,:N30_to_60_median]] = mapcols(zscore, md[:,[:S,:bSi,:TOC,:TIC,:yrBP,:N30_to_60_median]])
#md[!,Not([:Label,:period,:d])] = mapcols(f, md[:,Not([:Label,:period,:d])])

metadata_path = "./NetworkAnalysis/data/metadata.csv"
CSV.write(metadata_path, md)


# Step 4: run flashweave and save results
    netw = learn_network(data_path, metadata_path, max_k=1, n_obs_min=10, sensitive=true, heterogeneous=false)
    save_network("./NetworkAnalysis/output/network_$tag.jld2", netw)
# easy access
netw = load_network("./NetworkAnalysis/output/network_$tag.jld2")

    ##### Additionally, construct three networks with k=0  with sp in at least 3 and 5 samples ########
    k=1;#k=0
    netw2 = learn_network(data_path, metadata_path, max_k=k, n_obs_min=10, sensitive=true, heterogeneous=false)
    netw3 = learn_network(data_path3, metadata_path, max_k=k, n_obs_min=10, sensitive=true, heterogeneous=false)
    netw5 = learn_network(data_path5, metadata_path, max_k=k, n_obs_min=10, sensitive=true, heterogeneous=false)
    save_network("./NetworkAnalysis/output/network_k$(k)_min2_$tag.jld2", netw2)
    save_network("./NetworkAnalysis/output/network_k$(k)_min3_$tag.jld2", netw3)
    save_network("./NetworkAnalysis/output/network_k$(k)_min5_$tag.jld2", netw5)
    ##########################################################################################



