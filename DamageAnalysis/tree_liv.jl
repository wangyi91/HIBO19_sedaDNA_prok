#!/usr/bin/env julia
using DataFrames, JLD2, Pipe
using FlashWeave
using Graphs, SimpleWeightedGraphs
using Colors
using Compose, GraphPlot

include("./MicrobeProfiling/_get_annotated_df.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2")
netw = load_network("./NetworkAnalysis/output/network_$tag.jld2")
g = graph(netw)
cliper=clique_percolation(g,k=3)
ntaxa = sum(netw.meta_variable_mask.==0)

# load liv, ht from ../DmgMixtureModel/living_species.jl

# make taxonomic trees
using DelimitedFiles
writedlm("./MicrobeProfiling/output/taxlist_liv_Bacteria.csv", String.(unique(filter(:tax_path=>p->occursin("Bacteria",p),liv).tax_name)))
writedlm("./MicrobeProfiling/output/taxlist_liv_Archaea.csv", String.(unique(filter(:tax_path=>p->occursin("Archaea",p),liv).tax_name)))

# send taxlist_xxx.csv file to local
rsync -avP tvg137@dandycomp01fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/taxlist_liv*.csv /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/MicrobeProfiling/output/

# use phyloT to make trees, save in newick format and move to
mv /Users/yiwang/Downloads/liv_Bacteria*.newick /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Bacteria/
mv /Users/yiwang/Downloads/liv_Archaea*.newick /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Archaea/

# use R to remove bootstrap values and collapse internal nodes
# /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/tree_dominant_species.R

# send output trees back to julia
rsync -avP /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/*a/liv_*.modified.newick tvg137@dandycomp03fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output

using Phylo
ftax="Bacteria"; clrscheme=:Dark2_8; rank=" Phyla";
ftax="Archaea"; clrscheme=:tab10; rank=" Classes";
file="/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/liv_$(ftax).modified.newick"
tree = open(parsenewick, Phylo.path(file));
leafnames = getleafnames(tree);
#ph = PhyloBranches(tree)

df = @pipe DataFrame(nodename=getnodenames(tree)) |> 
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"p__"), n, missing))=>:phylum) |>
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"c__"), n, missing))=>:class)

lab = ftax[1]=="B" ? filter(!ismissing, df.phylum) : filter(!ismissing, df.class)
trait=map_depthfirst((val, node) -> val + ifelse(node in lab, findfirst(==(node),lab), 0), 0, tree, Int64);




# load hclust cutree result
group = DataFrame(taxname=unique(liv.tax_name), nodegroup=cutree(ht, k=4))
group = @pipe leftjoin(group, select(liv, [:tax_name,:tax_path])|>unique, on=:taxname => :tax_name) |>
         transform(_,:tax_path=>ByRow(p->tax_at_rank(p,"phylum"))=>:phylum) 
@pipe combine(groupby(group,[:phylum,:nodegroup]), nrow=>:n) |>sort(_,[:nodegroup,:n])



# now annotate the tree with key community groupings
df.nodegroup .=0;
for i in 1:maximum(group.nodegroup)
    splist = @pipe "s__".* filter(:nodegroup=>n->n==i, group).taxname|> replace.(_," "=>"_")
    df[findall(in(splist),df.nodename),:nodegroup] .= i
end

using Plots
colors = Dict(1=>cgrad(:tab10,categorical=true)[1], # blue
              2=>cgrad(:tab10,categorical=true)[3], # green
              3=>cgrad(:tab10,categorical=true)[2], # orange
              4=>cgrad(:tab10,categorical=true)[4], # red
              0=>:white)
alphas= Dict([i=>1 for i in 1:maximum(group.nodegroup)]); alphas[0]=0;
getindex.(Ref(colors), df.nodegroup)



# plot tree
gr(margins = 4Plots.cm, size=(2000, 20*length(leafnames)+200));
Plots.plot(tree, treetype=:dendrogram, showtips=false, legend=false,
           line_z = trait, linecolor=cgrad(clrscheme,rev=true), linewidth=0.5*ifelse(ftax=="Bacteria",1,2),
           markercolor=getindex.(Ref(colors), df.nodegroup), markeralpha=getindex.(Ref(alphas), df.nodegroup), markershape=:circle, markersize=4, markerstrokewidth=0);
savefig("./MicrobeProfiling/output/tree_liv_$(ftax).pdf");

# plot tree legends
using CairoMakie
fig=Figure(resolution=(1500,1000));
hidedecorations!(Axis(fig[1,1]));
hidedecorations!(Axis(fig[2,1]))
group_color = [MarkerElement(marker=:circle, color = clr, strokecolor = :transparent) for clr in collect(values(colors))];
Legend(fig[1,1], group_color, string.(collect(keys(colors))), "Groups (tree leaves)", tellheight = false,
              patchsize=(13,13), labelsize=13)

lrank = ftax[1]=="B" ? filter(!ismissing,df.phylum) : filter(!ismissing,df.class)
phylaclrs = palette(clrscheme, 1:length(lrank), rev=true)
phyla_color = [PolyElement(color = clr, strokecolor = :transparent, points = Point2f[(0,0.2),(1,0.2),(1,0.6),(0,0.6)]) for clr in phylaclrs];
Legend(fig[2,1], phyla_color, lrank, ftax*rank*" (tree branches)", nbanks=2, patchsize=(13,13), labelsize=13)
save("./MicrobeProfiling/output/tree_liv_$(ftax)_legends.pdf",fig)














