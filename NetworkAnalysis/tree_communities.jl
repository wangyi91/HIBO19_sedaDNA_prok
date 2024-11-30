#!/usr/bin/env julia
using DataFrames, JLD2, Pipe
using FlashWeave
using Graphs, SimpleWeightedGraphs
using Colors
using Compose, GraphPlot
#import Cairo, Fontconfig

#include("./MicrobeProfiling/_get_annotated_df.jl")

#pip="amaw"; alg=".lca";add="_ANI92";
#tag=pip*alg*add

#otudata = load_object("./InitialExploration/data/$tag.jld2")
#netw = load_network("./NetworkAnalysis/output/network_$tag.jld2")
#g = graph(netw)
#cliper=clique_percolation(g,k=3)
#ntaxa = sum(netw.meta_variable_mask.==0)
#cliq_idx = findall(>(2), length.(cliper));

# run first: ../NetworkAnalysis/key_communities.jl

nds = vcat(collect.(vcat(collect.([cliper[Is[i]] for i in 1:length(Is)])...))...) |> unique # full tree
nds = vcat(collect.(vcat(collect.([cliper[Is[i]] for i in [1,2,5,6,7]])...))...) |> unique # selected tree

splist = collect(skipmissing(tax_at_rank.(names(netw)[nds], "species")))
sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _)

# make taxonomic trees
using DelimitedFiles
writedlm("./NetworkAnalysis/output/taxlist_comm_Bacteria.csv", String.(unique(filter(:tax_path=>p->occursin("Bacteria",p),sub).tax_name)))
writedlm("./NetworkAnalysis/output/taxlist_comm_Archaea.csv", String.(unique(filter(:tax_path=>p->occursin("Archaea",p),sub).tax_name)))

# send taxlist_xxx.csv file to local
rsync -avP tvg137@dandycomp01fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/NetworkAnalysis/output/taxlist_comm*.csv /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/NetworkAnalysis/output/

# use phyloT to make trees, save in newick format and move to
mv /Users/yiwang/Downloads/comm_Bacteria*.newick /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Bacteria/
mv /Users/yiwang/Downloads/comm_Archaea*.newick /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Archaea/
# mv <path to tree>/comm_Bacteria*.newick <path to R proj>/R.pd/data/Bacteria/


# use R to remove bootstrap values and collapse internal nodes
# /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/tree_dominant_species.R

# send output trees back to julia
rsync -avP /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Bacteria/comm_*.modified.newick tvg137@dandycomp02fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/NetworkAnalysis/output
rsync -avP /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Archaea/comm_*.modified.newick tvg137@dandycomp02fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/NetworkAnalysis/output

using Phylo
ftax="Bacteria"; clrscheme=:tab20; rank="phylum";# :Dark2_8
ftax="Archaea"; clrscheme=:Accent_8; rank="phylum";
file="/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/NetworkAnalysis/output/comm_$(ftax).full.modified.newick"
tree = open(parsenewick, Phylo.path(file));
leafnames = getleafnames(tree);
#ph = PhyloBranches(tree)

df = @pipe DataFrame(nodename=getnodenames(tree)) |> 
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"p__"), n, missing))=>:phylum) |>
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"c__"), n, missing))=>:class)

phyla = filter(!ismissing, df.phylum);
class=filter(!ismissing, df.class);

ftax=="Bacteria" ? lab=phyla : lab=class; 
trait=map_depthfirst((val, node) -> val + ifelse(node in lab, findfirst(==(node),lab), 0), 0, tree, Int64);


# now annotate the tree with key community groupings
df.nodegroup .="";
df.ndgr .=0;
for i in [1,2,5,6,7]
#for i in 1:length(Is)
    clqs = collect.(cliper[cliq_idx])[Is[i]];
    clq = reduce(vcat, clqs) |> sort |> unique;
    splist = @pipe "s__".*collect(skipmissing(tax_at_rank.(names(netw)[clq], "species"))) |> replace.(_," "=>"_");
    df[findall(in(splist),df.nodename),:nodegroup] .= titles[i];
    df[findall(in(splist),df.nodename),:ndgr] .= i;
end
tmp=@pipe select(otudata,[:tax_name, :tax_path, :reference])|> unique |>
          transform(_,:tax_path=>ByRow(p->"s__".*replace(tax_at_rank(p,"species")," "=>"_"))=>:nodename)
anno = @pipe filter(:nodegroup=>g->g!="",df) |> select(_,[1,5]) |>
             leftjoin(_,tmp,on=:nodename) |> 
             transform(_,:tax_path=>ByRow(p->tax_at_rank(p,rank))=>Symbol(rank)) |> select(_,[1,2,5,6])

using Plots, Colors
rankpalette = DataFrame(Symbol(rank) => unique(anno[:,rank]), :rankcolor => "#".*hex.(palette(clrscheme, 1:length(unique(anno[:,rank])))))
grouppalette = DataFrame(ndgr=Vector(1:length(Is)), groupcolor=["#2986cc","#2CA02B","#ffdf4d","#744700","#901fb4","#FD7F0E","#D82728"]) 
#blue,green,yellow, brown,purple,orange,red
anno = @pipe leftjoin(anno, rankpalette, on=Symbol(rank)) |> leftjoin(_,grouppalette, on=:ndgr)

# here generate anno file for iTOL
anno[:,"Symbol"].=1; anno[:,"Size"].=1; anno[:,"Fill"].=1; anno[:,"Position"].=1;
rename!(anno, "nodename"=>"#nodename")

# for group/node annotation
CSV.write("./NetworkAnalysis/output/group_anno_$ftax.full.csv",anno[!,[1,7,8,6,9,10,2]])
# for rank/strip annotation
CSV.write("./NetworkAnalysis/output/rank_anno_$ftax.full.csv",anno[!,[1,5,4]])
# output accession if needed
CSV.write("./NetworkAnalysis/output/accession_$ftax.full.csv",anno[!,[1,3]])

# use bash to create annotation files for iTOL
cat ./data/iTOL_dataset_symbols.txt ./output/group_anno_Bacteria.full.csv > ./output/iTOL_group_annotation_Bacteria.full.txt
cat ./data/iTOL_dataset_strip.txt ./output/rank_anno_Bacteria.full.csv > ./output/iTOL_rank_annotation_Bacteria.full.txt

cat ./data/iTOL_dataset_symbols.txt ./output/group_anno_Archaea.full.csv > ./output/iTOL_group_annotation_Archaea.full.txt
cat ./data/iTOL_dataset_strip.txt ./output/rank_anno_Archaea.full.csv > ./output/iTOL_rank_annotation_Archaea.full.txt



















#=
CSV.write("./MicrobeProfiling/output/comm_groups_$ftax.csv",anno)
# now upload tree file (.newick) and anno file (.csv) to iTOL to generate annotated tree


colors = Dict(titles[1]=>cgrad(:tab10,categorical=true)[1], # blue
              titles[2]=>cgrad(:tab10,categorical=true)[3], # green
              titles[3]=>cgrad(:tab10,categorical=true)[5], # purple
              titles[4]=>cgrad(:tab10,categorical=true)[2], # orange
              titles[5]=>cgrad(:tab10,categorical=true)[4], # red
              titles[6]=>cgrad(:tab10,categorical=true)[6], # brown
              ""=>:white, "rest"=>:white)
alphas= Dict([titles[i]=>1 for i in 1:length(Is)]); alphas[""]=0;alphas["rest"]=0;
getindex.(Ref(colors), df.nodegroup)



# plot tree
gr(margins = 4Plots.cm, size=(2000, 20*length(leafnames)+200));
Plots.plot(tree, treetype=:fan, showtips=false, legend=false,
           line_z = trait, linecolor=cgrad(clrscheme,rev=true), linewidth=0.5*ifelse(ftax=="Bacteria",1,2),
           markercolor=getindex.(Ref(colors), df.nodegroup), markeralpha=getindex.(Ref(alphas), df.nodegroup), markershape=:circle, markersize=2, markerstrokewidth=0);
savefig("./MicrobeProfiling/output/tree_comm_$(ftax).pdf");

# plot tree legends
using CairoMakie
fig=Figure(resolution=(1500,1000));
hidedecorations!(Axis(fig[1,1]));
hidedecorations!(Axis(fig[2,1]))
group_color = [MarkerElement(marker=:circle, color = clr, strokecolor = :transparent) for clr in collect(values(colors))];
Legend(fig[1,1], group_color, collect(keys(colors)), "Groups (tree leaves)", tellheight = false,
              patchsize=(13,13), labelsize=13)

phyla = filter(!ismissing,df.phylum);
phylaclrs = palette(clrscheme, 1:length(phyla), rev=true)
phyla_color = [PolyElement(color = clr, strokecolor = :transparent, points = Point2f[(0,0.2),(1,0.2),(1,0.6),(0,0.6)]) for clr in phylaclrs];
Legend(fig[2,1], phyla_color, phyla, ftax*rank*" (tree branches)", nbanks=2, patchsize=(13,13), labelsize=13)
save("./MicrobeProfiling/output/tree_comm_$(ftax)_legends.pdf",fig)














