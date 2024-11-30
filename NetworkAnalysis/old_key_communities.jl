#!/usr/bin/env julia
using DataFrames, JLD2, Pipe
using FlashWeave, Graphs, SimpleWeightedGraphs
using Compose, GraphPlot, Colors

### This script outputs plot and table of low-abundance communities  ###

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2")
netw = load_network("./NetworkAnalysis/output/network_k1_min2_$tag.jld2")
g = graph(netw)
cliper=clique_percolation(g,k=3)
ntaxa = sum(netw.meta_variable_mask.==0)

    ## New grouping method: first based on connected components, then the rest based on depths
    titles = ["Deposited by floods","Early-Holocene assemblages","Present since Middle Ages","Present since modern time","New in recent years"]
    Is=[[],[],[],[],[]]; npanel = length(Is); # key community-groups; each element of Is is a vector of i's; i's are indices of cliq_idx
    cliq_idx = findall(>(2), length.(cliper));
    for (i,clq) in enumerate(collect.(cliper[cliq_idx]))
        splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
        sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _)
        if maximum(sub.Middle_depth) <= 10.5 push!(Is[5],i) 
        #if issubset([0.0,10.5], sub.Middle_depth) && length(unique(sub.Middle_depth))<=3 push!(Is[5],i) 
        elseif length(filter(d->d<122, unique(sub.Middle_depth)))>1 && length(filter(d->d>=122, unique(sub.Middle_depth))) <2 push!(Is[4],i)
        elseif length(filter(d->d<=300.5, unique(sub.Middle_depth)))>1 && length(filter(d->d>300.5, unique(sub.Middle_depth))) <2 push!(Is[3],i)
        elseif sum([966.1, 900.5, 624.3, 350.5] .∈ [sub.Middle_depth]) >1 && length(unique(sub.Middle_depth))<=5 push!(Is[1],i)
        elseif length(filter(d->d>=605.5, unique(sub.Middle_depth)))>1 && length(filter(d->d<605.5, unique(sub.Middle_depth)))<2 push!(Is[2],i)
        end
    end
    # count n of species in each community-group
    nsps = [@pipe vcat(collect.(cliper[cliq_idx])[Is[i]]...) |> unique |> count(x->x<=ntaxa, _) for i in 1:length(Is)]


# Part 1: make damage-depth histogram plot, stacked by phylo dist histogram plot. Groups as panels.
phylodist = load_object("./NetworkAnalysis/output/phylodist_$(tag).jld2"); # needs update everytime netw is updated
depths=otudata.Middle_depth |> sort |> unique;

using CairoMakie, Colors
fig = Figure(resolution=(220*(npanel+0.8)+100,1000));
axs = [Axis(fig[2,i], limits=(0,0.33,-depths[end]-10,-depths[1]+200), yticks=-depths, xticklabelsize=13, yticklabelsize=13,
            xlabel="DNA damage", ylabel="depth (cm)", 
            xlabelpadding=0, ylabelpadding=0) for i in 1:npanel]; # on axs dmg-depth is plotted, plus α and phylogenetic diversity
hidedecorations!.(axs, label=false, ticklabels = false, grid=false);
hideydecorations!.(axs[2:end], grid = false);
#hidexdecorations!.(axs[2:end], label=false, ticklabels=false, grid=false);
linkyaxes!(axs...);

bxs = [Axis(fig[1,i], limits=(0,4,nothing,nothing), xlabel="phylogenetic distance", ylabel="no. connections", xticklabelsize=13, yticklabelsize=13,
            xlabelpadding=0, ylabelpadding=0, title=titles[i]*"\n"*string(nsps[i])*" species") for i in 1:npanel]; # on bxs hist of phylo distance is plotted
hidedecorations!.(bxs, label=false, ticklabels = false);
rowsize!(fig.layout, 2, Relative(4/5));
colgap!(fig.layout, 10); rowgap!(fig.layout, 10);

cxs = [Axis(fig[2,i], limits=(-100,230,-depths[end]-10,-depths[1]+200), yticks=-depths, 
            xlabel="species richness", ylabel="", xaxisposition=:top, xlabelsize=12,
            xlabelpadding=0, ylabelpadding=0) for i in 1:npanel]; # on cxs species richness are plotted
hideydecorations!.(cxs); hidexdecorations!.(cxs, ticklabels=true);

for i in 1:length(Is)                                               
    clqs = collect.(cliper[cliq_idx])[Is[i]]    
    clq = reduce(vcat, clqs) |> sort |> unique                      
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _)
    # plot dmg histogram by depth
    for dep in depths
        dt = @pipe filter(:Middle_depth => d->d==dep, sub); if nrow(dt)==0 continue end
        #w = Float64.(sqrt.(dt.N_reads))*0.6
        CairoMakie.hist!(axs[i], Float64.(dt.damage); bins=0:0.005:0.4, normalization=:none, offset=-dep,weights=fill(5,size(dt)[1]), color=(:teal,0.5))
    end
    # plot species richness
    summ = combine(groupby(sub,:Middle_depth),nrow=>"n_species")
    deps = @pipe filter(d->d>=minimum(sub.Middle_depth), depths) |> filter(d->d<=maximum(sub.Middle_depth), _)
    app = DataFrame(Middle_depth=filter(d->!(d in sub.Middle_depth),deps)); app.n_species.=0; 
    @pipe append!(summ,app) |> sort!(_,:Middle_depth)
    CairoMakie.lines!(cxs[i], summ.n_species, summ.Middle_depth*(-1); linewidth=3, color=(:tan1,0.6))
    # plot phylo dist
    deleteat!(clq, clq .> ntaxa);# remove MVs
    wsub = Graphs.weights(g)[clq,clq]; w=vec(wsub);
    psub = phylodist[clq,clq]
    a=vec(psub); a=a[findall(>(0),w)];a=a[findall(in(1e-10..5),a)]
    #b=vec(psub); b=b[findall(<(0),w)];b=b[findall(in(1e-10..5),b)]
    CairoMakie.hist!(bxs[i], a; bins=0:0.1:5, normalization=:none, color=(:firebrick, 0.6))
    #CairoMakie.hist!(bxs[i], b; bins=0:0.1:5, normalization=:none, color=(:royalblue3, 0.3))
end
connect_color = [PolyElement(color = (clr, 0.4), strokecolor = :transparent) for clr in [colorant"firebrick", colorant"royalblue3"]];
#Legend(fig[1,npanel+1], connect_color, ["positive","negative"], "Connection type",
#       tellheight = false, tellwidth = false,labelsize=13)
#colsize!(fig.layout, npanel+1, Relative(0.7/(npanel+0.7)));
bar_color=PolyElement(color = (:teal, 0.5), strokecolor = :transparent, points = Point2f[(0.25,0),(0.55,0),(0.55,1),(0.25,1)]);
#line_color=PolyElement(color = (:tan1, 0.6), strokecolor = :transparent, points = Point2f[(0,0.25),(1,0.25),(1,0.55),(0,0.55)]);
#Legend(fig[2,npanel+1], [bar_color], ["1000?"],"Read count",
#       tellheight = false, labelsize=13)
#colsize!(fig.layout, npanel+1, Relative(0.7/(npanel+0.7)));
save("./NetworkAnalysis/output/hist_dmg_groups_manual.pdf", fig); 


# Part 2: output a table showing species metadata in each group





















