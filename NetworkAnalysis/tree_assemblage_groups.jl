#!/usr/bin/env julia
using DataFrames, JLD2, Pipe
using FlashWeave
using Graphs, SimpleWeightedGraphs
using Colors
using Compose, GraphPlot

# run first: ../NetworkAnalysis/key_communities.jl

include("env.jl") # WORK_DIR

nds = vcat(collect.(vcat(collect.([cliper[Is[i]] for i in 1:length(Is)])...))...) |> unique # full tree
nds = vcat(collect.(vcat(collect.([cliper[Is[i]] for i in [1,2,5,6,7]])...))...) |> unique # selected tree

splist = collect(skipmissing(tax_at_rank.(names(netw)[nds], "species")))
sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _)

# make taxonomic trees
using DelimitedFiles
writedlm("./NetworkAnalysis/output/taxlist_comm_Bacteria.csv", String.(unique(filter(:tax_path=>p->occursin("Bacteria",p),sub).tax_name)))
writedlm("./NetworkAnalysis/output/taxlist_comm_Archaea.csv", String.(unique(filter(:tax_path=>p->occursin("Archaea",p),sub).tax_name)))

# copy taxlist_xxx.csv file to local
`rsync -avP ./NetworkAnalysis/output/taxlist_comm*.csv <local_path>/NetworkAnalysis/output/`

# use phyloT to make trees, save in newick format and move to
`mv <path to tree>/comm_Bacteria*.newick <path to R proj>/R.pd/data/Bacteria/`
`mv <path to tree>/comm_Archaea*.newick <path to R proj>/R.pd/data/Archaea/`


# use ./R.pd/tree_assemblage_groups.R to remove bootstrap values and collapse internal nodes

# cp modified trees back 
`rsync -avP ./R.pd/data/Bacteria/comm_*.modified.newick ./NetworkAnalysis/output`
`rsync -avP ./R.pd/data/Archaea/comm_*.modified.newick ./NetworkAnalysis/output`

using Phylo
ftax="Bacteria"; clrscheme=:tab20; rank="phylum";# :Dark2_8
ftax="Archaea"; clrscheme=:Accent_8; rank="phylum";
file=WORK_DIR*"NetworkAnalysis/output/comm_$(ftax).modified.newick"
tree = open(parsenewick, Phylo.path(file));
leafnames = getleafnames(tree);

df = @pipe DataFrame(nodename=getnodenames(tree)) |> 
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"p__"), n, missing))=>:phylum) |>
   DataFrames.transform(_,:nodename=>ByRow(n->ifelse(startswith(n,"c__"), n, missing))=>:class)

phyla = filter(!ismissing, df.phylum);
class=filter(!ismissing, df.class);

ftax=="Bacteria" ? lab=phyla : lab=class; 
trait=map_depthfirst((val, node) -> val + ifelse(node in lab, findfirst(==(node),lab), 0), 0, tree, Int64);


# now annotate the tree with key assemblage groupings
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


