#!/usr/bin/env julia

# This script is to plot mapping statistics for taxa of interest 
using DataFrames
using CSVFiles, CSV
using StatsPlots
using Pipe: @pipe

# run ../NetworkAnalysis/key_communities.jl so that objects needed are loaded

push!(titles,"") # (update titles when ngroups is odd number)
push!(nsps, 0) # update number of species

tb=otudata


# Choose one of the four for plotting
## dmg vs ani
y="read_ani_mean";
y_range = 92:0.05:100 # ani_range
ylab="ANI %"; tit = "ani";

## dmg vs bratio
y="breadth_exp_ratio"
y_range = 0.73:0.001:1.01 # bratio_range
ylab="Breadth ratio"; tit="bratio";

## dmg vs cov_evenness
y="cov_evenness"
y_range = 0:0.001:0.813 # cov_evenness_range
ylab="Coverage evenness"; tit="covevenness";

## dmg vs coverage_covered_mean
y="coverage_covered_mean";
y_range = 1:0.01:5.3 # coverage_covered_mean_range
ylab="Mean covered depth"; tit="covcovedmean";



using KernelDensity, CairoMakie

dmg_range = 0:0.005:0.4
dens = kde((Vector{Float64}(tb.damage), Vector{Float64}(tb[!,y])));
bg = pdf(dens,dmg_range,y_range);

# Plot setup
fig = Figure(size=(1000,600))
axs = [Axis(fig[m,n], limits=(0,0.3,minimum(y_range),maximum(y_range)), title=permutedims(reshape(titles, (4,2)))[m,n]*"\n$(permutedims(reshape(nsps,(4,2)))[m,n]) taxa") for m in 1:2, n in 1:4]
ylabel = Label(fig[1:2, 0], ylab, rotation = pi/2)
xlabel = Label(fig[3,1:4], "DNA damage")

# Plot background and scatter data
for m in 1:2, n in 1:4
    i=m*(m+1)+n-2; if i==8 break end
    clqs = collect.(cliper[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe tb |> filter(:tax_name => n-> n in splist, _)
    # subset mapping statistics to taxa of interest
    tab = @pipe tb |> filter(:reference=>r->r in unique(sub.reference), _)
    # plot background (all taxa)
    CairoMakie.contourf!(axs[m,n], dmg_range, y_range, bg, levels = 10, colormap=:Greys_4)
    # plot taxa in the assemblage group
    CairoMakie.scatter!(axs[m,n], Vector{Float64}(tab.damage), tab[!,y],
                        color=tab.Middle_depth./maximum(tb.Middle_depth), colormap=cgrad(:starrynight, rev=true), colorrange = (0, 1),
                        markersize=5, alpha=0.7, strokewidth=0.1, strokecolor=(:black,0.3))
end


# Save the figure
save("./MicrobeProfiling/output/dmg_vs_$(tit)_assemblages_no_legends.pdf", fig)

# Create a separate figure for the colorbars
colorbar_fig = Figure(size=(200, 200)) # Adjust size as needed

# Add the contour plot colorbar
colorbar_contour = Colorbar(
    colorbar_fig[1, 1],
    colormap=:Greys_4,
    label="Frequency density \nof all taxa in\nthe dataset (%)",
    colorrange=(0, 100),
    vertical=true,
    height=70,
    labelrotation=0 
)

# Add the scatter plot colorbar
colorbar_scatter = Colorbar(
    colorbar_fig[2, 1],
    colormap=cgrad(:starrynight, rev=true),
    label= "Sample depth\n(cm)",
    colorrange=(0, 1283.8),
    vertical=true,
    height=70,
    labelrotation=0 
)

# Save the colorbar figure
save("./MicrobeProfiling/output/legends_for_dmg_vs_$(tit).pdf", colorbar_fig)



