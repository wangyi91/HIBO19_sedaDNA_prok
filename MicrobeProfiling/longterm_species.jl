#!/usr/bin/env julia

using DataFrames, StatsBase
using Pipe: @pipe

# This script is to analyse long-persisting species in four steps
include("./MicrobeProfiling/_get_annotated_df.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2")


# Step 1: visualise species abundance distribution
using CairoMakie
fig=Figure(resolution=(1300,450))
axs1=Axis(fig[1,1], xlabel="Species relative read abundance per sample",ylabel="Probability", xscale=log10, yticks=0:0.2:1)
axs2=Axis(fig[1,2], xlabel="Number of samples a species occurs in",ylabel="Species count")
axs3=Axis(fig[1,3], xlabel="Number of samples a species occurs in",ylabel="Species read count", yscale=log10)

#CairoMakie.hist!(axs1, otudata.N_reads, bins=1:2e3:6e4)
s=combine(groupby(otudata, :Label), :N_reads=>sum=>:sum_reads)
tmp=@pipe leftjoin(otudata, s, on=:Label) |> select(_, [:Label, :tax_id,:tax_name,:N_reads,:sum_reads]) |> 
    transform(_, [:N_reads,:sum_reads]=>ByRow((n,s)->n/s)=>:p_reads)
    CairoMakie.lines!(axs1,sort(tmp.p_reads), (1:nrow(tmp))./nrow(tmp))

tmp2 = @pipe otudata |> combine(groupby(_,:tax_name), nrow=>:noccur);
CairoMakie.hist!(axs2, tmp2.noccur, bins=1:1:34)

tmp3 = @pipe leftjoin(otudata, tmp2, on=:tax_name);
CairoMakie.boxplot!(axs3, tmp3.noccur, Int.(tmp3.N_reads))
save("./MicrobeProfiling/output/plot_SAD_rev.pdf",fig)




# Step 2 define long-persisting species and plot
keep = @pipe otudata|> combine(groupby(_,:tax_path), nrow=>:N_occur) |> filter(:N_occur=>o->o>20, _);
lp = @pipe otudata |> filter(:tax_path=>p->p in keep.tax_path);
save_object("./MicrobeProfiling/output/lp.jld2", lp)

# plot abundance vs depth heatmap with R.lda/heatmap_lp.R
wide = @pipe lp |> filter(:tax_path=>p->occursin(tax,p),_) |> unstack(_, :Middle_depth, :tax_name, :N_reads, fill=0)
CSV.write("./MicrobeProfiling/output/lp.csv",wide)
#mv /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/MicrobeProfiling/output/lp*.csv  /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.lda/Data/

# plot damage vs depth heatmap with R.lda/heatmap_lp.R
depths = unique(lp.Middle_depth);
fig = Figure(resolution=(400,1000))
Axis(fig[1,1], limits=(0,0.4,-depths[end]-10,-depths[1]+60), yticks=-depths, xticklabelsize=13, yticklabelsize=13,
            xlabel="DNA damage", ylabel="depth (cm)", xlabelpadding=0, ylabelpadding=0)
for dep in depths
        dt = @pipe filter(:Middle_depth => d->d==dep, lp); if nrow(dt)==0 continue end
        CairoMakie.hist!(Float64.(dt.damage); bins=0:0.005:0.4, normalization=:none, offset=-dep,
                         weights=fill(0.5,size(dt)[1]), color = (:teal, 0.5))
end
save("./MicrobeProfiling/output/hist_dmg_lp.pdf", fig);



# Step 3 plot overrepresented families, bac and arc separated
tax="Bacteria"
#tax="Archaea"
summ_lp=@pipe filter(:tax_path=>p->occursin(tax,p),lp) |> 
        select(_,[:tax_path,:tax_name,:read_ani_mean,:N_reads]) |>
        combine(groupby(_,[:tax_path,:tax_name]), :read_ani_mean=>minimum=>:min_ani,:N_reads=>mean=>:mean_reads) |> 
        transform(_,:tax_path=>ByRow(p->tax_at_rank(p,"phylum"))=>:phylum) |>
        transform(_,:tax_path=>ByRow(p->tax_at_rank(p,"family"))=>:family)

CSV.write("./MicrobeProfiling/output/table_lp_$tax.csv",rename(summ_lp,:mean_reads=>:avg_reads_per_occurrence))

a=@pipe combine(groupby(summ_lp,[:family, :phylum]), nrow=>:nsp,:mean_reads=>sum=>:nreads)
b=@pipe otudata|>select(_,:tax_path)|>unique|> transform(_,:tax_path=>ByRow(p->tax_at_rank(p,"family"))=>:family) |>combine(groupby(_,:family), nrow=>:nsp)
c=@pipe leftjoin(a,b,on=:family,makeunique=true)|>transform(_,[:nsp,:nsp_1]=>ByRow((a,b)->a/b)=>:freq) |>sort(_,:nreads) |>                 
                     filter([:nsp,:nreads]=>(s,r)->s>1 && r>3e3,_)# s>0 && r>50 for archaea
d=@pipe combine(groupby(c,:phylum),nrow=>:nfamily,:nsp=>sum=>:nsp,:nsp_1=>sum=>:nsp_1,:nreads=>sum=>:nreads)|> transform(_,[:nsp,:nsp_1]=>ByRow((a,b)->a/b)=>:freq) |>sort(_,:nreads,rev=true)

fig=Figure(resolution=(800, nrow(c)*15+100))
Axis(fig[1,1], yticks = collect(1:nrow(c)), ytickformat=tks-> string.([c.family[trunc(Int,tk)] for tk in tks]),
    xgridvisible=false,ygridvisible=false)
ylims!(low = 0,high=nrow(c)+1);xlims!(low = 0)
barplot!(collect(1:nrow(c)),c.nreads,direction=:x,flip_labels_at=2e4,color_over_bar=:white,# 1e4 for arc
         bar_labels=string.(c.nsp).*" (" .* string.(round.(c.freq.*100,digits=1)) .*"%)")
save("./MicrobeProfiling/output/barplot_lp_$(tax[1:3]).pdf",fig)


# a few filtering criteria
ksp1 = @pipe lp|>combine(groupby(_,:tax_path),nrow=>:noccur)|>filter(:noccur=>n->n>26,_).tax_path
ksp = @pipe filter(:Middle_depth=>d->d<=84.5,lp)|>combine(groupby(_,:tax_path),nrow=>:noccur)|>filter(:noccur=>n->n<1,_).tax_path
ksp = @pipe combine(groupby(lp,:tax_path),:Middle_depth=>minimum=>:min_depth)|>filter(:min_depth=>n->n==34.5,_).tax_path
# summarise and display
tmp = @pipe filter(:tax_path=>p->p in intersect(ksp,ksp1),lp) |> combine(groupby(_,[:tax_path,:tax_name]),:Middle_depth=>minimum=>:min_depth,:N_reads=>mean=>:mean_reads,:read_ani_mean=>minimum=>:min_ani)|> transform(_,:tax_path=>ByRow(p->tax_at_rank(p,"family"))=>:family)|>transform(_,:tax_path=>ByRow(p->tax_at_rank(p,"phylum"))=>:phylum)|>sort(_,[:min_depth,:phylum,:family])|> filter(:mean_reads=>n->n>300,_)

# output above df as tables


# Step 4 plot mapping stats with ./plot_mapping_stats_lp.jl

