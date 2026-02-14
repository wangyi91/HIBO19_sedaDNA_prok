#!/usr/bin/env julia
using DataFrames, JLD2, Pipe
using FlashWeave, Graphs, SimpleWeightedGraphs
using Compose, GraphPlot, Colors
using CairoMakie, Colors

### This script outputs plot and table of low-abundance communities  ###

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")
include("./NetworkAnalysis/_spectral_clustering.jl")
include("./NetworkAnalysis/_plot_communities.jl") # contains function generate_plot()
include("./NetworkAnalysis/_heatmap_table_of_groups.jl") # contains function output_tb_hmp_of_groups()

pip="amaw"; alg=".lca";add="_ANI92";rank="species";
tag=pip*alg*add

otudata = load_object("./deContamination/data/$tag.$rank.jld2")
netw = load_network("./NetworkAnalysis/output/network_$tag.jld2")
g = graph(netw)
cliper=clique_percolation(g,k=3)
cliq_idx = findall(>(2), length.(cliper));
ntaxa = sum(netw.meta_variable_mask.==0)

phylodist = load_object("./NetworkAnalysis/output/phylodist_$(tag).jld2"); # needs update everytime netw is updated

## Grouping of clique-percolation generated communities: spectral clustering
ngroups=7 
grouping = sp_cluster(cliper, cliq_idx, otudata, ngroups)
ntax = [@pipe vcat(collect.(cliper[cliq_idx])[grouping[i]]...) |> unique |> count(x->x<=ntaxa, _) for i in 1:length(grouping)]

is=grouping[sortperm(ntax, rev=true)[1:7]] # choose the 7 most taxa abundant groups to present
display_order = [6,5,7,1,3,4,2]; # reorder Is for data presentation
Is=is[display_order]

titles = ["Long-term, stable",
          "Long-term, unstable",
          "Holocene uncommon 1",
          "Holocene uncommon 2",
          "Linked to flood depositions",
          "Predominant since the Middle Ages", 
          "Predominant since modern time"]

# count n of species in each community-group
nsps = [@pipe vcat(collect.(cliper[cliq_idx])[Is[i]]...) |> unique |> count(x->x<=ntaxa, _) for i in 1:length(Is)]

# Plot: make damage-depth histogram plot, stacked by phylo dist histogram plot. Groups as panels.
generate_plot()

# Output source data table
output_hmp_tb_of_groups()

