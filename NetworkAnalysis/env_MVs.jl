#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe
using FlashWeave
using Graphs, SimpleWeightedGraphs

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2");k=0;#k=1
netw2 = load_network("./NetworkAnalysis/output/network_k$(k)_min2_$tag.jld2");
netw3 = load_network("./NetworkAnalysis/output/network_k$(k)_min3_$tag.jld2");
netw5 = load_network("./NetworkAnalysis/output/network_k$(k)_min5_$tag.jld2");

####################### Part 1: Most connected MVs ##############################################################
dfs = [DataFrame(),DataFrame(),DataFrame()]

for (i,netw) in enumerate([netw2, netw3,netw5])
    g=graph(netw)
    # degree of nodes, only positive connections
    pc = sum(Graphs.weights(g).>0, dims=2) |> vec;
    # degree of nodes, only negative connections
    nc = sum(Graphs.weights(g).<0, dims=2) |> vec;
    # a summary dataframe of node connections
    cnodes = DataFrame(node_name=names(netw), degree=degree(g), pos_connect=pc, neg_connect=nc);
    # find the most connected continous MVs
    #mvs = findall(x->netw.meta_variable_mask[x]==1 && !startswith(names(netw)[x], "d_")
    #              && !any(occursin.(["bSi","TOC","TIC","S","N30_to_60_median","yrBP"],names(netw)[x])),
    #              1:length(netw.variable_ids));
    mvs = findall(x->netw.meta_variable_mask[x]==1 && !startswith(names(netw)[x], "d_"), 1:length(netw.variable_ids));

    dfs[i] = @pipe cnodes[mvs,1:4] |> sort(_,:degree);
    dfs[i].noccur .= maximum([2i-1,2])
end

#join the three tables    
df = @pipe append!(dfs[1],dfs[2]) |> append!(_,dfs[3]);

# plot
using DataStructures
landtypes = SortedDict("11"=>"Urban",
                 "12"=>"Mixed settlements",
                 "23"=>"Rainfed villages",
                 "24"=>"Pastoral villages",
                 "32"=>"Residential rainfed croplands",
                 "33"=>"Populated croplands",
                 "41"=>"Residential rangelands",
                 "42"=>"Populated rangelands",
                 "43"=>"Remote rangelands",
                 "51"=>"Residential woodlands",
                 "52"=>"Populated woodlands",
                 "53"=>"Remote woodlands",
                 "54"=>"Inhabited drylands",
                 "61"=>"Wild woodlands",
                 "62"=>"Wild drylands");
geodata = SortedDict("S"=>"sulfur",
                    "TIC"=>"total inorganic carbon",
                   "TOC"=>"total organic carbon",
                  "bSi"=>"biogenic silica",
                 "N30_to_60_median"=>"temperature",
                 "yrBP"=>"continuous time");
climdata = SortedDict("N30_to_60_median"=>"temperature","yrBP"=>"time (continous)");
periods = SortedDict("period_bronze_age"=>"Bronze Age",
 "period_contemporary"=>"Contemporary",
 "period_early_modern"=>"Early Modern",
 "period_iron_age"=>"Iron Age",
 "period_late_modern"=>"Late Modern",
 "period_mesolithic"=>"Mesolithic",
 "period_middle_ages"=>"Middle Ages",
 "period_neolithic"=>"Neolithic",
 "period_roman"=>"Roman Time")

dicts=merge(landtypes,geodata,climdata,periods);


tit=["Land Use","landuse"]
#tit=["(Pre)historic Periods","historic"]
tit=["Environmental Variables","envs"]
if tit[1]=="Land Use" 
    summ=filter(:node_name=>n->!startswith(n,"p"),df) 
elseif tit[1]=="Land Use"
    summ =filter(:node_name=>n->startswith(n,"p"),df) 
else summ=filter(:node_name=>n->n in(["bSi","TOC","TIC","S","N30_to_60_median","yrBP"]),df)
end
#add mv_id for plotting
ticks=summ.node_name |> unique;
DataFrames.transform!(summ, :node_name=>ByRow(n->findall(==(n), ticks)[1])=>:mv_id);
summ=stack(summ,3:4);
summ.ctype.=(summ.variable.=="pos_connect")*1;

ytl=unique(summ[!,[:node_name,:mv_id]]); # a df for yticklabels

# order by chronological order instead of no. connections
k==0 ? ytl.mv_id=[3,7,4,2,5,1,6,8,9] : ytl.mv_id=[7,9,4,1,6,2,3,5,8] 
summ=leftjoin(select(summ,Not(:mv_id)),ytl, on=:node_name)


using CairoMakie
fig=Figure(resolution=(1000,nrow(ytl)*50+50))
axs = Axis(fig[1, 1], title=tit[1], xlabel="Number of connections", ylabel="", yticks=ytl.mv_id,
           ytickformat=tks->[dicts[ytl[ytl.mv_id.==tk,:].node_name[1]] for tk in tks]);

CairoMakie.barplot!(axs, summ.mv_id, summ.value, direction=:x,stack=summ.ctype,
                    dodge=ifelse.(summ.noccur.<5, 5 .- summ.noccur, mod.(6,summ.noccur)),
                    color=cgrad(:tab20;categorical=true, rev=true)[ifelse.(summ.noccur.==2,summ.noccur.-1,summ.noccur).+summ.ctype])

shades = [:grey50,:grey90];
colors = [cgrad(:tab20;categorical=true, rev=true)[i] for i in [2,4,6]];
connect_shade = [PolyElement(strokecolor = :transparent, color=sh) for sh in shades];
netw_color = [PolyElement(color = clr, strokecolor = :transparent) for clr in colors];
Legend(fig[1,2], [connect_shade, netw_color], [["positive","negative"],string.([2,3,5])], ["Type of connection", "Species minimum occurences\nin network"])

save("./NetworkAnalysis/output/barplot_MV_k$(k)_$(tit[2]).pdf", fig);





##################  NOT USING  ##### Part 2: Most connected species ##############################################
# first rerun Part 1 first 4 commands to get cnodes for netw1
ntaxa = sum(netw.meta_variable_mask.==0)

# Calculate inter- and intra-rank connections for each species
# calculate if connections are within a given rank
rank = "phylum" # rank to evaluate
# construct a matrix showing if connections are intra-rank (1) or inter-rank (0)
nm = names(netw)[1:ntaxa]
rnm = tax_at_rank.(nm,rank)
inrank = broadcast(==, rnm, permutedims(rnm))

# add to df "cnodes"
cnodes.pos_connect_intra_phylum.=0;
cnodes.neg_connect_intra_phylum.=0;
cnodes.pos_connect_intra_phylum[1:ntaxa] = @pipe (Graphs.weights(g).>0)[1:ntaxa,1:ntaxa] .* inrank |> sum(_;dims=2) |> vec;
cnodes.neg_connect_intra_phylum[1:ntaxa] = @pipe (Graphs.weights(g).<0)[1:ntaxa,1:ntaxa] .* inrank |> sum(_;dims=2) |> vec;

cnodes.pos_connect_inter_phylum.=0;
cnodes.neg_connect_inter_phylum.=0;
cnodes.pos_connect_inter_phylum[1:ntaxa] = @pipe (Graphs.weights(g).>0)[1:ntaxa,1:ntaxa] .* (inrank.==0) |> sum(_;dims=2) |> vec;
cnodes.neg_connect_inter_phylum[1:ntaxa] = @pipe (Graphs.weights(g).<0)[1:ntaxa,1:ntaxa] .* (inrank.==0) |> sum(_;dims=2) |> vec;

# keep only species
cnodes = @pipe filter(:node_name=>n->occursin("s__",n),cnodes)|> sort(_,:degree)

# plot species degree distribution
fig=Figure()
axs = Axis(fig[1, 1], xlabel="Degree of species", xticks=0:8:48)
CairoMakie.hist!(axs, cnodes.degree, bins=0:4:52)
save("./NetworkAnalysis/output/degree_distribution_species.pdf", fig);




