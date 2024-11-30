#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe
using CSV, CSVFiles
using Diversity, Phylo
using FlashWeave
using Graphs, SimpleWeightedGraphs
using Compose, GraphPlot, Colors
import Cairo, Fontconfig

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

netw = load_network("./NetworkAnalysis/output/network_$tag.jld2")
g = graph(netw)

# Part 1: phylogenetic distance between taxa of interest
#
############## Run only when tree is updated using tree_all.jl #################
using Phylo
file1 = "/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/Bacteria_all_all_50.modified.newick"
tree1 = open(parsenewick, Phylo.path(file1))
#leafnames = getleafnames(tree1)
nodenames1 = getnodenames(tree1)

file2 = "/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/Archaea_all_all_50.modified.newick"
tree2 = open(parsenewick, Phylo.path(file2))
nodenames2 = getnodenames(tree2)


ntaxa = sum(netw.meta_variable_mask.==0)
netw_nodenames = names(netw)[1:ntaxa]

nodenames_std = "s__".*replace.(tax_at_rank.(netw_nodenames,"species"), " "=>"_")

# calculate phylogenetic distance between two taxa
phylodist = zeros(ntaxa, ntaxa)
Threads.@threads for j in 1:ntaxa-1
    for i in j+1:ntaxa
        if in(nodenames_std[i],nodenames1)&&in(nodenames_std[j], nodenames1)
            phylodist[i,j] = distance(tree1, nodenames_std[i], nodenames_std[j])
        elseif in(nodenames_std[i],nodenames2)&&in(nodenames_std[j], nodenames2)
            phylodist[i,j] = distance(tree2, nodenames_std[i], nodenames_std[j])
        else phylodist[i,j] = 999
        end
    end
end

save_object("./NetworkAnalysis/output/phylodist_$tag.jld2", phylodist)
#############################################################

phylodist = load_object("./NetworkAnalysis/output/phylodist_$tag.jld2")

a=vec(phylodist)

adjmt = adjacency_matrix(g)[1:ntaxa,1:ntaxa]
# pos
sub = (adjmt.>0) .* phylodist 
b=vec(sub)
w = (adjmt.>0).* adjmt .* (phylodist .>0) 

using StatsPlots
gr(margins = 1.5Plots.cm)
StatsPlots.density(b[findall(in(1e-10..10),b)], weights = collect(w[findall(in(1e-10..10),w)]), labels="positive");
StatsPlots.density!(a[findall(in(1e-10..10),a)], labels="random");
savefig("./NetworkAnalysis/output/hist_phylodist_positive$add3.pdf");

# neg
sub = (adjmt.<0) .* phylodist 
b=vec(sub)
w = (adjmt.<0).* adjmt .* (phylodist .>0) * (-1)

gr(margins = 1.5Plots.cm)
StatsPlots.density(b[findall(in(1e-10..10),b)], weights = collect(w[findall(in(1e-10..10),w)]), labels="negative");
StatsPlots.density!(a[findall(in(1e-10..10),a)], labels="random");
savefig("./NetworkAnalysis/output/hist_phylodist_negative$add3.pdf");






