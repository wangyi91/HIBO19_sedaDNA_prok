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
gr(margins = 1.5Plots.cm, size=(850, 500))
StatsPlots.boxplot(df.community_size, df.N_reads, fill=:skyblue,yaxis=:log10, marker=(1,:black),line=(1,:black),
               xlabel="Community size", ylabel="Read count of each species occurence", label="nreads", legend = false);
StatsPlots.scatter!(twinx(),summ.community_size, summ.n_tax, yaxis=:log10, markercolor=:orange, 
                    ylabel="Total number of unique species", label="ntax", legend = false);
savefig("./NetworkAnalysis/output/boxplot_scomm_vs_nreads.pdf");





gr(margins = 1.5Plots.cm, size=(1000, 500))
StatsPlots.boxplot(df.community_size, df.N_reads, fill=:skyblue,yaxis=:log10, marker=(1,:black),line=(1,:black), 
               xlabel="community size", ylabel="read count of each species occurence", label="nreads", legend=:outertopright);
StatsPlots.scatter!(twinx(),summ.community_size, summ.n_tax, yaxis=:log10, markercolor=:orange, ylabel="number of unique species", label="ntax", legend=:outerright);
#@df summ plot!(twinx(),:community_size, :n_tax,yaxis=:log10
#               ylabel="number of unique species", linecolor=:orange, label="");
#@df summ StatsPlots.scatter!(twinx(),:community_size, :n_tax, markercolor=:orange, ylabel="number of unique species", label="ntax", legend=:outerright);
savefig("./NetworkAnalysis/output/boxplot_scomm_vs_nreads_local.pdf");


# Plot 2: number and tax composition of singletons in each sample
single = @pipe single |> DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank(p, "phylum"))=> :phylum) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "class"))=> :class) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "family"))=> :family) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "genus"))=> :genus)

summ = @pipe combine(groupby(single, [:Label, :Middle_depth, :phylum]), nrow=>"n_species") |>
             DataFrames.transform(_, :phylum=>ByRow(r->findall(==(r), unique(_.phylum))[1])=>:stack)
ticks=summ.Middle_depth |> unique |> reverse
summ = @pipe summ |> DataFrames.transform(_, :Middle_depth=>ByRow(d-> findall(==(d), ticks)[1])=>:depth_id)

using CairoMakie
fig=Figure(resolution=(1250,1000))
colsize!(fig.layout, 1, Relative(2/3))
axs = Axis(fig[1, 1], yticks=unique(summ.depth_id), ytickformat=tks-> string.([summ[summ.depth_id.==tk,:Middle_depth][1] for tk in tks]), 
          ylabel="Depth (cm)", xlabel="Number of singleton species")
hidedecorations!(axs, label=false, ticklabels = false);

CairoMakie.barplot!(axs, summ.depth_id, summ.n_species, colormap=:seaborn_muted, color=summ.stack, stack=summ.stack, direction=:x)

colors = cgrad(:seaborn_muted, [0,1])[unique(summ.stack)./maximum(summ.stack)];
group_color = [PolyElement(color = color, strokecolor = :transparent) for color in colors];
Legend(fig[1,2], group_color, unique(summ.phylum), "phylum", nbanks = 2)
save("./NetworkAnalysis/output/barplot_singleton.pdf", fig);


# Plot 3: taxonomic tree of singleton species, coloured by the depth they appear
# make taxonomic trees
    using DelimitedFiles
    writedlm("./MicrobeProfiling/output/taxlist_single100_Bacteria.csv", 
             String.(unique(filter([:tax_path,:N_reads]=>(p,n)->occursin("Bacteria",p)&&n>99,single).tax_name)))
    writedlm("./MicrobeProfiling/output/taxlist_single50_Archaea.csv", String.(unique(filter(:tax_path=>p->occursin("Archaea",p),single).tax_name)))

    # send taxlist_xxx.csv file to local
    rsync -avP tvg137@dandycomp01fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/taxlist_single*.csv /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/MicrobeProfiling/output/

    # use phyloT to make trees, save in newick format and move to R.pd
    mv /Users/yiwang/Downloads/single100_Bacteria.newick /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Bacteria/
    mv /Users/yiwang/Downloads/single50_Archaea.newick /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Archaea/

    # use R to remove bootstrap values and collapse internal nodes
    # /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/tree_dominant_species.R

    # send output trees back to julia
    rsync -avP /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Bacteria/single100_Bacteria.modified.newick tvg137@dandycomp01fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output

using Phylo
ftax="Bacteria"; clrscheme=:Dark2_8; rank=" Phyla";
ftax="Archaea"; clrscheme=:tab10; rank=" Classes";
file="/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/single100_$(ftax).modified.newick"
tree = open(parsenewick, Phylo.path(file));
leafnames = getleafnames(tree);

df = @pipe DataFrame(nodename=getnodenames(tree)) |>
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"p__"), n, missing))=>:phylum) |>
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"c__"), n, missing))=>:class)

phyla = filter(!ismissing, df.phylum);
class=filter(!ismissing, df.class);

ftax=="Bacteria" ? lab=phyla : lab=class;
trait=map_depthfirst((val, node) -> val + ifelse(node in lab, findfirst(==(node),lab), 0), 0, tree, Int64);


# now annotate the tree with depth
df.nodegroup .=9999.9;
for d in unique(single.Middle_depth)
    splist = @pipe filter(:Middle_depth=>dep->dep==d, single).tax_name |> replace.(_," "=>"_")
    df[findall(in(splist),df.nodename),:nodegroup] .= d
end

using Plots
colors = Dict([d=>cgrad(:matter, [0,1])[d/maximum(unique(single.Middle_depth))] for d in unique(single.Middle_depth)]);
colors[9999.9] = RGBA{Float64}(1.0,1.0,1.0,1.0);# white

alphas= Dict([d=>1 for d in unique(single.Middle_depth)]); alphas[9999.9]=0;
getindex.(Ref(colors), df.nodegroup)

# plot tree
gr(margins = 4Plots.cm, size=(2000, 20*length(leafnames)+200));
Plots.plot(tree, treetype=:fan, showtips=false, legend=false,
           line_z = trait, linecolor=cgrad(clrscheme,rev=true), linewidth=0.5*ifelse(ftax=="Bacteria",1,2),
           markercolor=getindex.(Ref(colors), df.nodegroup), markeralpha=getindex.(Ref(alphas), df.nodegroup), markershape=:circle, markersize=5, markerstrokewidth=0);
savefig("./MicrobeProfiling/output/tree_single_$(ftax).pdf");



























# plot for each individual cliq community (larger than a certain size)
cliq_idx = findall(>(10), length.(cliper));

ntaxa = sum(netw.meta_variable_mask.==0)
using CairoMakie#Plots

for (i,clq) in enumerate(collect.(cliper[cliq_idx]))
    len = length(clq)
    # plot graph of clip communities (highlighted) with neighbours
 #=   member = vcat(taxmember.idx,fill(4,sum(netw.meta_variable_mask)));  
    node_hl = clq; member[node_hl].=5;
    println(i);println(node_hl);
    node_hl_neb = vcat(collect.([all_neighbors(g, n) for n in node_hl])...) |> unique |> sort;
    sg, vmap = induced_subgraph(g, node_hl_neb)

    member_sg = member[vmap];
    ndfill_sg = nodecolor[member_sg];
    ndsize_sg = [1,1,1,2,1][member_sg];
    ndlabel_sg = hide_labels_w_str(names(netw)[vmap], ["__"]);
    wts_sg = [weights(sg)[src.(edges(sg))[i], dst.(edges(sg))[i]] for i in 1:length(edges(sg))];
    strcolor_sg = strokecolor[(wts_sg.>0).+1];
=#
    #sgp = gplot(sg, nodefillc=ndfill_sg, nodesize=ndsize_sg, nodelabel=ndlabel_sg, edgestrokec=strcolor_sg, layout=spring_layout);
    #draw(PDF("./NetworkAnalysis/output/subgraph$(add3)_cliq_neb_$(i)_size$len.pdf", 16cm, 16cm), sgp);

    splist = collect(names(netw)[clq]) #splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], ftax)))
    sub = @pipe otudata |> filter(:tax_path => n-> n in splist, _)
    
    # plot histogram of pairwise phylo distances for a clique community
    #deleteat!(clq, clq .> ntaxa);# remove MVs
    #psub = phylodist[clq,clq]
    #gr(margins = 1.5Plots.cm)
    #a=vec(psub); a=a[findall(>(0),a)];a=a[findall(<(100),a)]
    #histogram(a, bins=0:0.1:5)
    #savefig("./NetworkAnalysis/output/phylodist_cliq_$(i)_size$len.pdf");

    # plot density of dmg at each depth
    fig = Figure(resolution=(500,1000))
    depths=otudata.Middle_depth |> sort |> unique
    Axis(fig[1, 1], title = "damage distribution",limits=(0,0.4,-depths[end]-100,-depths[1]+200), yticks=-depths,
              xlabel="damage", ylabel="depth (cm)", yreversed=false)
    for dep in depths
        dt = @pipe filter(:Middle_depth => d->d==dep, sub); if nrow(dt)==0 continue end
        w = Float64.(sqrt.(dt.N_reads))
        d = dt.Middle_depth[1]
        CairoMakie.density!(Float64.(dt.damage), npoints = 5000, weights = w, offset=-d,
                            color = (:lightseagreen, 0.4), bandwidth = 0.002)
    end
    save("./NetworkAnalysis/output/hist_dmg_$(tag)_$(i)_size$(len).pdf", fig);
end





# plot histogram of N_reads for the largest clique communties
i=1;
for clq in collect.(clique_percolation(g, k=4)[cliq_idx])
    splist = "s__".*collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _)
    gr(margins = 1.5Plots.cm)
    histogram(sqrt.(sub.N_reads))
    savefig("./NetworkAnalysis/output/cliq_Nreads_$i.pdf")

i+=1
end

