#!/usr/bin/env julia

using DataFrames
using Pipe: @pipe
using JLD2

# This script is to plot mapping statistics for long-persisting and living species
pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add

otudata = load_object("./InitialExploration/data/$tag.jld2")

# load lp (from ./MicrobeProfiling/long_persisting_species.jl) and liv (from ./DmgMixtureModel/living_species.jl)
lp=load_object("./MicrobeProfiling/output/lp.jld2")
liv=load_object("./DmgMixtureModel/output/liv.jld2")
nsps=length.(unique.([lp.tax_path,liv.tax_path]))

using KernelDensity
dmg_range = 0:0.005:0.4
# Choose one of the four for plotting
## dmg vs ani
y="read_ani_mean";
y_range = 92:0.05:100 # ani_range
ylab="ANI %"; tit = "ani";

## dmg vs bratio
y="breadth_exp_ratio"
y_range = 0.75:0.001:1 # bratio_range
ylab="Breadth ratio"; tit="bratio";

## dmg vs cov_evenness
y="cov_evenness"
y_range = 0:0.001:0.813 # cov_evenness
ylab="Coverage evenness"; tit="covevenness";

## dmg vs coverage_covered_mean
y="coverage_covered_mean";
y_range = 1:0.01:5.3 # coverage_covered_mean_range
ylab="Mean covered depth"; tit="covcovedmean";

dens = kde((Vector{Float64}(otudata.damage), Vector{Float64}(otudata[!,y])));
bg = pdf(dens,dmg_range,y_range);


# set up plot
using CairoMakie
fig = Figure(resolution=(800,300))
axs = [Axis(fig[1,i], limits=(nothing,nothing,minimum(y_range),maximum(y_range)), title=["Long-term species\n$(nsps[1]) species","Sediment-dwelling species\n$(nsps[2]) species",""][i]) for i in 1:3]
ylabel = Label(fig[:, 0], ylab, rotation = pi/2, tellheight=false)
xlabel = Label(fig[2,1:2], "DNA damage")

for i in 1:2
    i==1 ? tb=lp : tb=liv
    # plot background (all taxa)
#    CairoMakie.contourf!(axs[i], dmg_range, y_range, bg, levels = 10, colormap=:Greys_4)
    # plot species in community-group
    CairoMakie.scatter!(axs[i], Vector{Float64}(tb.damage), tb[!,y], 
                        color=tb.Middle_depth./maximum(otudata.Middle_depth), colormap=cgrad(:starrynight, rev=true), colorrange = (0,1), 
                        markersize=4, alpha=0.5, strokewidth=0.001, strokecolor=(:black,0.3))
end
#Colorbar(fig[1,4],sc)
#Legend(fig[1,4],[],[],["Depth","Density"],tellheight = false,)
save("./MicrobeProfiling/output/dmg_vs_$(tit)_lpliv.pdf", fig);


