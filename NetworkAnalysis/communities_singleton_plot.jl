#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe
using FlashWeave
using Graphs, SimpleWeightedGraphs
using Compose, GraphPlot, Colors
#import Cairo, Fontconfig

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2")
netw = load_network("./NetworkAnalysis/output/network_$tag.jld2")
g = graph(netw)
cliper=clique_percolation(g,k=3)

# Plot 1: cliq community size vs N_reads and n_tax
com=DataFrame();
for clq in collect.(cliper) #clique_percolation(g, k=3)
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    tmp = @pipe otudata |> filter(:tax_name=>n-> n in splist,_)
    tmp.community_size .= length(clq)
    com = vcat(tmp, com, cols=:union)
end
# singleton
singleton = @pipe otudata |> combine(groupby(_,:tax_name), nrow=>"N_occur") |> filter(:N_occur=>o->o==1, _);
single = @pipe otudata |> filter(:tax_name=>n-> n in singleton.tax_name,_);
single.community_size .= -4;
# not in k-clique communities
not = @pipe otudata |> filter(:tax_name => n-> !(n in singleton.tax_name) && !(n in com.tax_name),_);
not.community_size .= 0;
df = vcat(com,not,single, cols=:union);


summ = @pipe df|> unique(_,[:tax_name,:community_size]) |> combine(groupby(_, :community_size), nrow=>:n_tax)

using StatsPlots
default(fontfamily = "Helvetica")
gr(margins = 1.5Plots.cm, size=(850, 500))
StatsPlots.boxplot(df.community_size, df.N_reads, fill=:skyblue,yaxis=:log10, marker=(1,:black),line=(1,:black),
               xlabel="Size of assemblages identified with clique percolation", ylabel="Read count of each taxon occurence", label="nreads", legend = false);
StatsPlots.scatter!(twinx(),summ.community_size, summ.n_tax, yaxis=:log10, markercolor=:orange, 
                    ylabel="Total number of unique taxa", label="ntax", legend = false);
savefig("./NetworkAnalysis/output/boxplot_scomm_vs_nreads.pdf");





