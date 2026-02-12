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
frank="species"


# Step 1: load data_training by running ./InitialExploration/vcat.jl 
otudata=load_object("./deContamination/data/$tag.$frank.jld2")

# Step 2: filter out singletons, convert data to otu table and save as csv
keep = @pipe otudata |> combine(groupby(_, :tax_path), nrow=>"N_occur") |> filter(:N_occur=>o->o>=2, _)

otu = @pipe otudata |> filter(:tax_path => p->p in keep.tax_path, _) |>
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


# Step 5: plot
g = graph(netw)
cliper=clique_percolation(g,k=3)
#save_object("./NetworkAnalysis/output/cliper_$tag.jld2",cliper)





taxmember = @pipe DataFrame(tax_path=names(otu)[2:end]) |>
        DataFrames.transform(_, :tax_path => ByRow(t -> tax_at_rank.(t, "kingdom")) => :kingdom) |>
        DataFrames.transform(_, :kingdom => ByRow(t -> ifelse(t=="Bacteria",1,ifelse(t=="Archaea",2,3))) => :idx)

member0 = vcat(taxmember.idx,fill(4,sum(netw.meta_variable_mask)));
member=copy(member0)
nodecolor = [colorant"turquoise",colorant"orange",colorant"lightgrey",colorant"cornsilk4",colorant"palevioletred"]; # bacteria, archaea, others, MVs, highlighted
    #nodecolor = [nothing,nothing,colorant"lightgrey",colorant"cornsilk4",colorant"palevioletred"]; # selective coloring
ndfill = nodecolor[member];
ndsize = [1,1,1,2,1][member];
ndlabel = hide_labels_w_str(names(netw), ["__",":\"","period"]);

strokecolor = [colorant"cornflowerblue",colorant"lightgrey"];
    #strokecolor = [nothing,colorant"lightgrey"];
wts = [Graphs.weights(g)[src.(edges(g))[i], dst.(edges(g))[i]] for i in 1:length(edges(g))];
strcolor = strokecolor[(wts.>0).+1];
    # nodes size proportional to their degree
    #nodesize = [Graphs.outdegree(g, v) for v in Graphs.vertices(g)]
    #alphas = nodesize/maximum(nodesize)
    #nodefills= [RGBA(nodecolor[member[i]],1) for i in 1:length(alphas)]

gp = gplot(g, nodefillc=ndfill, nodesize=ndsize, NODESIZE=0.01, 
           nodelabel=ndlabel, NODELABELSIZE=2, nodelabeldist=2.5, nodelabelangleoffset=π/4, 
           edgestrokec=strcolor, layout=spring_layout);
draw(PDF("./NetworkAnalysis/output/graph$tag.pdf", 16cm, 16cm), gp)


# find cliques and color the nodes differently
cliper=clique_percolation(g, k=3)
save_object("./NetworkAnalysis/output/cliper_k3$(add3).jld2",cliper);

cliq_idx = findall(>(6), length.(clique_percolation(g, k=3)));
#cliq_idx = findall(in(6..10), length.(clique_percolation(g, k=4)));
node_hl = vcat(collect.(clique_percolation(g, k=3)[cliq_idx])...) |> unique |> sort;
member = member0; member[node_hl].=5; # hightlight these nodes with #5 color
ndfill = nodecolor[member]; # update ndfill


# Step 6: community analysis

## option 1: plot subgraph of communities detected using clique percolation 
sg, vmap = induced_subgraph(g, node_hl)

## option 2: plot subgraph of above selected communities + their direct neighbours
node_hl_neb = vcat(collect.([all_neighbors(g, n) for n in node_hl])...) |> unique |> sort;
sg, vmap = induced_subgraph(g, node_hl_neb)

member_sg = member[vmap];
ndfill_sg = nodecolor[member_sg];
ndsize_sg = [1,1,1,2,1][member_sg];
ndlabel_sg = hide_labels_w_str(names(netw)[vmap], ["__"]);
wts_sg = [weights(sg)[src.(edges(sg))[i], dst.(edges(sg))[i]] for i in 1:length(edges(sg))];
strcolor_sg = strokecolor[(wts_sg.>0).+1];

sgp = gplot(sg, nodefillc=ndfill_sg, nodesize=ndsize_sg, nodelabel=ndlabel_sg, edgestrokec=strcolor_sg, layout=spring_layout, NODESIZE=0.01, NODELABELSIZE=2.5);
draw(PDF("./NetworkAnalysis/output/subgraph$(add3)_cliq.pdf", 16cm, 16cm), sgp);





























# get subgraph based on node connection
connected_components(g)

for id in findall(>(10),length.(connected_components(g))) # nodes of interest
    subnodes = connected_components(g)[id] # largest subsets: 1,2,4,209
    #subnodes = findall(>(4),degree(g)) # degree > 4

    sg, vmap = induced_subgraph(g, subnodes)
    member_sg = member[vmap];
    ndfill_sg = nodecolor[member_sg];
    ndsize_sg = [1,1,1,2,1][member_sg];
    ndlabel_sg = hide_labels_w_str(names(netw)[vmap], ["__"]);
    wts_sg = [weights(sg)[src.(edges(sg))[i], dst.(edges(sg))[i]] for i in 1:length(edges(sg))]
    strcolor_sg = strokecolor[(wts_sg.>0).+1]

    sgp = gplot(sg, nodefillc=ndfill_sg, nodesize=ndsize_sg, nodelabel=ndlabel_sg, edgestrokec=strcolor_sg, NODESIZE=0.01, NODELABELSIZE=2.5, layout=spring_layout);
    draw(PDF("./NetworkAnalysis/output/subgraph$add3$id.pdf", 16cm, 16cm), sgp);
end

    # get original aDNA data of taxa from subgraph
    #x = @pipe data_training |> filter(:tax_name => n->n in names(netw)[vmap],_) |> filter(:Middle_depth => d->d<408,_)
    # => then use ./DmgMixtureModel/plot_dmg_density_on_single_plot.jl to make plot



# access the graph
betweenness_centrality(g)
connected_components(g)
dominating_set(g, DegreeDominatingSet())
independent_set(g, DegreeIndependentSet())
global_clustering_coefficient(g)
adjacency_matrix(g)==weights(g)

# plot histogram of weights
a=vec(adjacency_matrix(g))
histogram(a[findall(!=(0),a)])

# plot degree distribution
using Plots
histogram(degree(g));
savefig("./NetworkAnalysis/output/degree_distribution$add3.pdf");

# check degree of MVs
netw.variable_ids[end-64:end-34]
degree(g)[end-64:end-34]





