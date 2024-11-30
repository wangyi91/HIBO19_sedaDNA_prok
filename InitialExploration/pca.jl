#!/usr/bin/env julia


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
