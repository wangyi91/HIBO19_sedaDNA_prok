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
include("env.jl") # contains WORK_DIR

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

netw = load_network("./NetworkAnalysis/output/network_$tag.jld2")
g = graph(netw)

# --- Calculate pair-wise phylogenetic distance between taxa ---
#
# Trees are built with phyloT and then modified using R.pd/tree.R
using Phylo
file1 = WORK_DIR*"/HIBO_shotgun/analysis/MicrobeProfiling/output/Bacteria_all_all_50.modified.newick"
tree1 = open(parsenewick, Phylo.path(file1))
nodenames1 = getnodenames(tree1)

file2 = WORK_DIR*"HIBO_shotgun/analysis/MicrobeProfiling/output/Archaea_all_all_50.modified.newick"
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

# to load phylodist
phylodist = load_object("./NetworkAnalysis/output/phylodist_$tag.jld2")

