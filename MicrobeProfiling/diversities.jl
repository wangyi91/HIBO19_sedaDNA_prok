#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe

include("./MicrobeProfiling/_get_annotated_df.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add
rank="species"

otudata = load_object("./deContamination/data/$tag.$rank.jld2")

# plot relative read abundance by phylum in each sample, bac and arc separate
dt = @pipe otudata |> DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank(p, "phylum"))=> :phylum) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "class"))=> :class) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "family"))=> :family) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "genus"))=> :genus)

using Diversity.Hill, CairoMakie

pd=load_object("./MicrobeProfiling/output/pd.jld2")
taxa=["Bacteria","Archaea"];
fig=Figure(size=(800,500))
axs = [Axis(fig[i,1], limits=(nothing,nothing,0,nothing), title=taxa[i], ylabel="Shannon diversity",ylabelcolor=:dodgerblue3, 
            xlabel="Years BP", xticklabelsize=10, xticks=unique(otudata.yrBP),xtickformat=tks->[string(round(Int,tk)) for tk in tks],
            xreversed=true, xticklabelrotation=3.1415926/2) for i in 1:2]
bxs = [Axis(fig[i,1], limits=(nothing,nothing,0,[50,8][i]),
            ylabel="Phylogenetic diversity",ylabelcolor=:red3, yaxisposition = :right, xreversed=true) for i in 1:2]
hidexdecorations!.(bxs)
for (i,tax) in enumerate(taxa)
    div = @pipe dt |> select(_,[:tax_path,:N_reads,:Label,:Middle_depth, :yrBP]) |>
    filter(:tax_path => p -> occursin(tax, p), _) |>
    filter(:N_reads => n -> n > [500, 100][i], _) |>
    groupby(_, :Label) |>
    transform(_, :N_reads => (n -> n ./ sum(n)) => :rel_read_abd) |>
    combine(groupby(_, [:Middle_depth, :yrBP]),
            :rel_read_abd => p -> hillnumber(p, 1).diversity[1]) |>
    rename(_,:rel_read_abd_function=>:shannon)
#    div = @pipe dt|>filter(:tax_path=>p->occursin(tax,p))|> filter(:N_reads=>n->n>[100,50][i],_) |>
#            combine(groupby(_,[:Middle_depth,:yrBP]),:N_reads=>n->hillnumber(n,1).diversity[1])|> 
#            rename(_,:N_reads_function=>:shannon)
    div.pd=pd[!,i+1]
    scatterlines!(axs[i],div.yrBP, div.shannon, color=:dodgerblue3)
    scatterlines!(bxs[i],div.yrBP, div.pd, linestyle=:dash, color=:red3,markersize=8)
end
save("./MicrobeProfiling/output/lineplot_diversities.pdf",fig)

