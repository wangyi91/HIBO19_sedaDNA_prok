#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe

include("./MicrobeProfiling/_get_annotated_df.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2")

# plot relative read abundance by phylum in each sample, bac and arc separate
dt = @pipe otudata |> DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank(p, "phylum"))=> :phylum) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "class"))=> :class) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "family"))=> :family) |>
DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "genus"))=> :genus)

using Diversity.Hill, CairoMakie

pd=load_object("./MicrobeProfiling/output/pd.jld2")
taxa=["Bacteria","Archaea"];
fig=Figure()
axs = [Axis(fig[i,1], limits=(nothing,nothing,0,nothing), title=taxa[i], ylabel="Shannon diversity",ylabelcolor=:dodgerblue3, 
            xlabel="Years BP", xticklabelsize=12, xticks=unique(otudata.yrBP),xtickformat=tks->[string(round(Int,tk)) for tk in tks],
            xreversed=true, xticklabelrotation=π/2) for i in 1:2]
bxs = [Axis(fig[i,1], limits=(nothing,nothing,0,[50,8][i]),
            ylabel="Phylogenetic diversity",ylabelcolor=:red3, yaxisposition = :right, xreversed=true) for i in 1:2]
hidexdecorations!.(bxs)
for (i,tax) in enumerate(taxa)
    div = @pipe dt|>filter(:tax_path=>p->occursin(tax,p))|> filter(:N_reads=>n->n>[100,50][i],_) |>
            combine(groupby(_,[:Middle_depth,:yrBP]),:N_reads=>n->hillnumber(n,1).diversity[1])|> 
            rename(_,:N_reads_function=>:shannon)
    div.pd=pd[!,i+1]
    scatterlines!(axs[i],div.yrBP, div.shannon, color=:dodgerblue3)
    scatterlines!(bxs[i],div.yrBP, div.pd, linestyle=:dash, color=:red3,markersize=8)
end
save("./MicrobeProfiling/output/lineplot_diversities.pdf",fig)






fig=Figure()
axs = Axis(fig[1,1])
lines!(Div.Middle_depth, Div.shannon)
save("./MicrobeProfiling/output/lineplot_diversity.pdf",fig)

using CairoMakie
fig=Figure(resolution=(1250,1000))
#colsize!(fig.layout, 1, Relative(2/3))
axs = Axis(fig[1, 1], yticks=unique(summ.depth_id), ytickformat=tks-> string.([summ[summ.depth_id.==tk,:Middle_depth][1] for tk in tks]), 
          ylabel="Depth (cm)", xlabel="Relative read counts")
hidedecorations!(axs, label=false, ticklabels = false);

CairoMakie.barplot!(axs, summ.depth_id, summ.preads, colormap=:seaborn_muted, color=summ.stack, stack=summ.stack, direction=:x)

colors = cgrad(:seaborn_muted, [0,1])[unique(summ.stack)./maximum(summ.stack)];
group_color = [PolyElement(color = color, strokecolor = :transparent) for color in colors];
Legend(fig[1,2], group_color, unique(summ.phylum), "phylum", nbanks = 2)
save("./MicrobeProfiling/output/barplot_diversity.pdf", fig);











tmp=combine(groupby(dt, :Label),:N_reads=>sum=>:Nreads)
summ = @pipe combine(groupby(dt, [:Label, :Middle_depth, :phylum]),:N_reads=>sum=>:nreads) |>
             leftjoin(_, tmp, on=:Label) |>
             transform(_, [:nreads,:Nreads]=>ByRow((n,N)->1+n/N)=>:preads) |>
             transform(_, :phylum=>ByRow(r->findall(==(r), unique(_.phylum))[1])=>:stack)
ticks=summ.Middle_depth |> unique |> reverse
summ = @pipe summ |> transform(_, :Middle_depth=>ByRow(d->findall(==(d),ticks)[1])=>:depth_id)

using CairoMakie
fig=Figure(resolution=(1250,1000))
colsize!(fig.layout, 1, Relative(2/3))
axs = Axis(fig[1, 1], yticks=unique(summ.depth_id), ytickformat=tks-> string.([summ[summ.depth_id.==tk,:Middle_depth][1] for tk in tks]), 
          ylabel="Depth (cm)", xlabel="Relative read counts")
hidedecorations!(axs, label=false, ticklabels = false);

CairoMakie.barplot!(axs, summ.depth_id, summ.preads, colormap=:seaborn_muted, color=summ.stack, stack=summ.stack, direction=:x)

colors = cgrad(:seaborn_muted, [0,1])[unique(summ.stack)./maximum(summ.stack)];
group_color = [PolyElement(color = color, strokecolor = :transparent) for color in colors];
Legend(fig[1,2], group_color, unique(summ.phylum), "phylum", nbanks = 2)
save("./MicrobeProfiling/output/barplot_diversity.pdf", fig);


