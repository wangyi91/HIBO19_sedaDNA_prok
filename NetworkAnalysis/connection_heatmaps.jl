#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe
using FlashWeave
using Graphs, SimpleWeightedGraphs
using Compose, GraphPlot, Colors

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

# This script uses objects from key_communities.jl 
pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2")
netw = load_network("./NetworkAnalysis/output/network_k1_min2_$tag.jld2")
g = graph(netw)
cliper=clique_percolation(g,k=3)
ntaxa = sum(netw.meta_variable_mask.==0)

# load Is, cliq_idx, nsps... from ./key_communities_plot.jl 

# Part 1: output number of connections for species in each community-group
using Clustering, Plots
for i in 1:length(Is)
    clqs = collect.(clique_percolation(g, k=3)[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    deleteat!(clq, clq .> ntaxa);
    wsub = Graphs.weights(g)[clq,clq]; # weight matrix of community
    pconn = sum(x->x>0, wsub, dims=1) |> vec; # n of +connections for each sp
    nconn = sum(x->x<0, wsub, dims=1) |> vec; # n of -connections for each sp
    tab = @pipe DataFrame(tax_path=names(netw)[clq],
                 kingdom=tax_at_rank.(names(netw)[clq],"kingdom"),
                 phylum=tax_at_rank.(names(netw)[clq],"phylum"),
                 class=tax_at_rank.(names(netw)[clq],"class"),
                 order=tax_at_rank.(names(netw)[clq],"order"),
                 family=tax_at_rank.(names(netw)[clq],"family"),
                 genus=tax_at_rank.(names(netw)[clq],"genus"),
                 species=tax_at_rank.(names(netw)[clq],"species"),
                 n_pos_conn = pconn, n_neg_conn = nconn, group=titles[i]) |> 
            leftjoin(_,select(otudata,[:tax_path,:Label,:Middle_depth,:N_reads]),on=:tax_path)|>sort(_,:n_pos_conn,rev=true)
    c=@pipe combine(groupby(tab,[:tax_path,:family,:order,:class,:phylum]),nrow=>:noccur, :N_reads=>mean=>:avg_reads)|>sort(_,[:noccur,:avg_reads], rev=true)        
    @pipe combine(groupby(c,[:order,:phylum]),nrow=>:n)|>sort(_,:n)
    CSV.write("./NetworkAnalysis/output/species_data_commgr$i.csv", tab)
    if i==1 Tab=tab else Tab=vcat(Tab,tab) end
end
CSV.write("./NetworkAnalysis/output/community_groups_species_data.csv", Tab)



# Part 2: plot connection (+ or -) heatmaps of each community-group at chosen taxonomic rank
# get list of rank names ordered as in taxonomic tree
using Phylo
ftax="Bacteria"; clrscheme=:Dark2_8;
#ftax="Archaea"; clrscheme=:Accent_8;
file="/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/$(ftax).modified.newick"
tree = open(parsenewick, Phylo.path(file));
leafnames = getleafnames(tree);

df = @pipe DataFrame(nodename=getnodenames(tree)) |>
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"s__"), n, missing))=>:species)

ords=filter(!ismissing, df.species);

# define a function which outputs the order of a list of rank names as in taxonomic tree
function ord_in_tree(sp; ords=ords)
    out=findfirst(==(replace.("s__"*sp," "=>"_")), ords)
    return(out)
end

ranks=["phylum","class","order","family","genus"];
rk=4; # so that ranks[rk] is plotted in heatmap

using Plots, DelimitedFiles
for i in 1:length(Is), ctype in ["+","-"]
    clqs = collect.(clique_percolation(g, k=3)[cliq_idx])[Is[i]]
    clq0 = reduce(vcat, clqs) |> sort |> unique
    deleteat!(clq0, clq0.>ntaxa);
    # prepare a df where species are sorted as in taxonomic tree [BACTERIA ONLY]
    sorted_clq = @pipe DataFrame(nd_id=clq0, tax_path=names(netw)[clq0], 
                     kingdom=tax_at_rank.(names(netw)[clq0],"kingdom"),
                     phylum=tax_at_rank.(names(netw)[clq0],"phylum"),
                     class=tax_at_rank.(names(netw)[clq0],"class"),
                     order=tax_at_rank.(names(netw)[clq0],"order"),
                     family=tax_at_rank.(names(netw)[clq0],"family"),
                     genus=tax_at_rank.(names(netw)[clq0],"genus"),
                     species=tax_at_rank.(names(netw)[clq0],"species")) |>
               filter(:kingdom=>d->d=="Bacteria",_) |> sort(_, order(:species, by=ord_in_tree))
    clq = sorted_clq.nd_id

    ctype == "+" ? w=(Graphs.weights(g).>0)[clq, clq] : w=(Graphs.weights(g).<0)[clq, clq]; # choose to output pos or neg connections
    mt = @pipe w |> Matrix |> DataFrame(_, :auto);
    mt.node_name = names(netw)[clq];
    DataFrames.transform!(mt, :node_name => ByRow(n->tax_at_rank(n, ranks[rk]))=>:node_rank);
    to_group = names(mt[!, Not([:node_rank,:node_name])]);
    mt = @pipe mt |> combine(groupby(_, :node_rank), to_group .=> sum .=> to_group) |> permutedims(_, 1);
    mt.node_name = names(netw)[clq];
    DataFrames.transform!(mt, :node_name => ByRow(n->tax_at_rank(n, ranks[rk]))=>:node_rank);
    to_group = names(mt[!, Not([:node_rank,:node_name])]);
    mt = @pipe mt |> combine(groupby(_, :node_rank), to_group .=> sum .=> to_group)

    if nrow(mt)>30 # avoid heatmap being too large
        cut=sort(maximum.(eachcol(mt[:,2:end])), rev=true)[30]
        rowkeep=findall(>(cut),maximum.(eachcol(mt[:,2:end])))
        mt=mt[rowkeep,append!([1],rowkeep.+1)]
    end
    using Plots    
    gr(margins = 2Plots.cm, size=(600, 600))
    plot(heatmap(mt.node_rank, mt.node_rank, Matrix(mt[:,2:end]), c=cgrad(ifelse(ctype=="+",:Reds, :Blues))),
         xticks = :all, xtickfontsize=5, xrotation = 90, yticks = :all, ytickfontsize=5, dpi=800);
    title!(titles[i]*"\n$(nsps[i]) species");
    savefig("./NetworkAnalysis/output/hmp_$(ctype)connect_$(ranks[rk])_commgr$i.png");
    # output a table annotated with necessary taxonomic ranks
    tab = @pipe leftjoin(mt, select(sorted_clq, Symbol.(ranks[1:rk]))|>unique, on=[:node_rank=>Symbol(ranks[rk])])|>
            select(_,append!([:node_rank],Symbol.(ranks[1:rk-1])))
    tab.id=1:nrow(tab); sort!(tab,:id,rev=true); select!(tab,Not(:id))
    CSV.write("./NetworkAnalysis/output/table_$(ctype)connect_$(ranks[rk])_commgr$i.csv", tab)
end
