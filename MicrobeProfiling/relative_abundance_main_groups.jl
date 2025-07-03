#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe

include("./MicrobeProfiling/_get_annotated_df.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add
frank="species"

otudata = load_object("./deContamination/data/$tag.$frank.jld2")

# plot relative read abundance by phylum in each sample, bac and arc separate
dt = @pipe otudata |> DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank(p, "phylum"))=> :phylum) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "class"))=> :class) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "family"))=> :family) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "genus"))=> :genus)


tmp=combine(groupby(dt, :Label),:N_reads=>sum=>:Nreads)

summ = @pipe dt |> filter(:Middle_depth => d-> d<20, _) |> 
             combine(groupby(_, [:Label, :Middle_depth, :class]),:N_reads=>sum=>:nreads) |>
             leftjoin(_, tmp, on=:Label) |>
             transform(_, [:nreads,:Nreads]=>ByRow((n,N) -> n/N)=>:relative_read_abd) |>
             transform(_, :class=>ByRow(r->findall(==(r), unique(_.class))[1])=>:stack)

@pipe summ |> filter(:relative_read_abd=>r->r>0.04,_) |> sort(_, :relative_read_abd) |> 
              filter(:class=>c->c in ["Phycisphaerae",
                                      "Actinobacteria",
                                      "Verrucomicrobiae",
                                      "Bacteroidia",
                                      "Verrucomicrobiae",
                                      "Anaerolineae"],_)
