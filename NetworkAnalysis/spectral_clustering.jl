#!/usr/bin/env julia
using DataFrames, JLD2, Pipe
using FlashWeave, Graphs, SimpleWeightedGraphs
using Compose, GraphPlot, Colors

### This script uses spectral clustering to group low-abundance communities  ###

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2")
netw = load_network("./NetworkAnalysis/output/network_k1_min2_$tag.jld2")
g = graph(netw)
cliper=clique_percolation(g,k=3)
ntaxa = sum(netw.meta_variable_mask.==0)

# prepare a table where rows are samples, columns are each clique (community), values are number of taxa
tab=DataFrame(Label=String[], sum_reads=Int64[], n_taxa=Int64[], N_clq=Int64[])
cliq_idx = findall(>(2), length.(cliper));
for (i,clq) in enumerate(collect.(cliper[cliq_idx]))
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _) |> 
                combine(groupby(_,:Label), :N_reads=>sum=>:sum_reads, nrow=>:n_taxa) |>
                insertcols(_,:N_clq=>i)
    append!(tab, sub)
end

wide = unstack(tab, :N_clq, :Label, :n_taxa, fill=0)

# Option 1: use spectral clustering to group communities
## output to local and use R
CSV.write("./NetworkAnalysis/output/matrix_communities.csv", wide)







using LinearAlgebra, Distances
  m = Matrix(wide[:,2:end]) # value = reads
  m = log1p.(wide[:,2:end]|> Matrix{Float64}) #value = log(reads)
  m = Matrix((wide[:,2:end] .>0) .* 1) # value = presence/absence (1/0)



# Option 1: use spectral clustering to group communities
## output to local and use R
#using DelimitedFiles
#writedlm("./NetworkAnalysis/output/matrix_communities.csv", m)
CSV.write("./NetworkAnalysis/output/matrix_communities.csv", wide)

## generate similarity matrix using Euc distance
simM = pairwise(euclidean, m, dims=1)
eigvals(simM)
eigvecs(simM)

# Option 2: use PCA to group communities (obsolete)  
using MultivariateStats
M = fit(PCA, m; maxoutdim=3)

pj = @pipe projection(M) |> DataFrame(_, :auto);
pj.N_clq=wide.N_clq

using Plots
h1=30;l1=40; h2=60;l2=40; h3=60;l3=20;
gr(margins = 0.5Plots.cm, size=(1000, 425));
p1=plot(pj.x1, pj.x2, pj.x3, marker_z=pj.N_clq,
     color=:roma, seriestype=:scatter, mode="markers",
     label="", camera = (h1,l1), legend=false, markerstrokewidth=0.1);
p2=plot(pj.x1, pj.x2, pj.x3, marker_z=pj.N_clq,
     color=:roma, seriestype=:scatter, mode="markers",
     label="", camera = (h2,l2), legend=false, markerstrokewidth=0.1);
p3=plot(pj.x1, pj.x2, pj.x3, marker_z=pj.N_clq,
     color=:roma, seriestype=:scatter, mode="markers",
     label="", camera = (h3,l3), legend=false,markerstrokewidth=0.1);
plot(p1,p2,p3,layout=@layout [p1 p2 p3]);
savefig("./NetworkAnalysis/output/pca_clqs.pdf");













### Run after loading data from local-lca.jl ###
# prepare a matrix
wide = unstack(dt, :yrBP, :tax_name, :N_reads, fill=0)
wide_norm = unstack(dt, :yrBP, :tax_name, :p_reads, fill=0)

# normalise by row (samples) to account for unequal depths for eukaryote reads
using LinearAlgebra
if pip=="holi"
    #m = mapslices(x -> x / norm(x), Matrix(wide[:,2:end]), dims=2)
    m = Matrix(wide_norm[:,2:end])
else
  # log transform, then to matrix
  m = log1p.(wide[:,2:end]|> Matrix{Float64})
end

using MultivariateStats
M = fit(PCA, m; maxoutdim=3)

pj = @pipe projection(M) |> DataFrame(_, :auto);
#pj.Middle_depth = wide.Middle_depth
pj.yrBP = wide_norm.yrBP;

using Plots
h1=30;l1=40; h2=60;l2=40; h3=60;l3=20;
gr(margins = 0.5Plots.cm, size=(1000, 425));
p1=plot(pj.x1, pj.x2, pj.x3, marker_z=pj.yrBP, 
     color=:roma, seriestype=:scatter, mode="markers",
     label="", camera = (h1,l1), legend=false, markerstrokewidth=0.1);
p2=plot(pj.x1, pj.x2, pj.x3, marker_z=pj.yrBP,
     color=:roma, seriestype=:scatter, mode="markers",
     label="", camera = (h2,l2), legend=false, markerstrokewidth=0.1);
p3=plot(pj.x1, pj.x2, pj.x3, marker_z=pj.yrBP,
     color=:roma, seriestype=:scatter, mode="markers",
     label="", camera = (h3,l3), legend=false,markerstrokewidth=0.1);
plot(p1,p2,p3,layout=@layout [p1 p2 p3]);
savefig("./InitialExploration/output/pca_$(pip)_$(tax)_$(frank).pdf");
