#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe
using FlashWeave
using Graphs, SimpleWeightedGraphs

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

dt = load_object("./InitialExploration/data/$tag.jld2");k=0;#k=1

netw2 = load_network("./NetworkAnalysis/output/network_k$(k)_min2_$tag.jld2");
netw3 = load_network("./NetworkAnalysis/output/network_k$(k)_min3_$tag.jld2");
netw5 = load_network("./NetworkAnalysis/output/network_k$(k)_min5_$tag.jld2");

####################### Part 1: Most connected MVs ##############################################################
dfs = [DataFrame(),DataFrame(),DataFrame()]

for (i,netw) in enumerate([netw2, netw3, netw5])
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
fig=Figure(size=(1000,nrow(ytl)*50+50))
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


####################### Part 2: Find taxa connected to MVs of interest ##############################################
nw = load_network("./NetworkAnalysis/output/network_k1_min5_$tag.jld2");
nw = load_network("./NetworkAnalysis/output/network_k0_min5_$tag.jld2");
g = graph(nw)

# Step 1: Find the vertex index for "N30_to_60_median"
v = findfirst(==( "N30_to_60_median" ), names(nw))

# Step 2: Get the indices of its neighbors
neighbor_indices = neighbors(g, v)

# Step 3: Map those indices back to names
neighbor_names = names(nw)[neighbor_indices]

neighbor_taxa = filter(!ismissing, tax_at_rank.(neighbor_names, "species"))

# Separate neighbors by positive/negative connection
pos_neighbors = Int[]
neg_neighbors = Int[]

for u in neighbor_indices
    w = Graphs.weights(g)[v, u]
    if w > 0
        push!(pos_neighbors, u)
    elseif w < 0
        push!(neg_neighbors, u)
    end
end

# Map indices to names
pos_neighbor_names = names(nw)[pos_neighbors]
neg_neighbor_names = names(nw)[neg_neighbors]
pos_neighbor_taxa = filter(!ismissing, tax_at_rank.(pos_neighbor_names, "species"))
neg_neighbor_taxa = filter(!ismissing, tax_at_rank.(neg_neighbor_names, "species"))

# now go to key_communities.jl and load clique percolation
all_matches = String[]  # will hold all collected taxa

for i in 1:length(Is)
    clqs = collect.(cliper[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    
    matches = intersect(pos_neighbor_taxa, splist)
    println(length(matches))
    append!(all_matches, matches)
end

# Drop missings
all_matches = filter(!ismissing, all_matches)

dt_sub = filter(:tax_name=>n->n in all_matches, dt)
sub = unstack(dt_sub,:Label, :tax_path, :N_reads)[:,2:end]


using CairoMakie
#sample_names = 1:nrow(sub)                # or actual sample IDs if you have them
sample_names = @pipe dt_sub.yrBP |> unique; sample_names = round.(Int, sample_names)
taxa_paths = names(sub)

n_taxa = length(taxa_names)
ncols = 5                            # how many plots across
nrows = ceil(Int, n_taxa / ncols)

fig = Figure(size = (ncols * 300, nrows * 300))

for (i, pth) in enumerate(taxa_paths)
    taxon = tax_at_rank.(pth,"species")
    counts = sub[:, pth] |> collect
    ax = Axis(fig[div(i-1, ncols)+1, mod(i-1, ncols)+1],
              title = "$(taxon)\n($(tax_at_rank(pth, "family")))", titlefont = :italic,
              xticks = (1:length(sample_names), string.(sample_names)), xticklabelsize = 8,
              xticklabelrotation = 3.1416/2, xticklabelalign = (:right, :center))
    
    barplot!(ax, 1:length(sample_names), counts, color = :steelblue)
end
save("./NetworkAnalysis/output/temp_taxa_barplot.png", fig)
