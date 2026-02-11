#!/usr/bin/env julia

using DataFrames
using Pipe: @pipe

# This script is to model ANI change through time for high-ANI long-persisting species
# It continues after first loading "otudata" from ./NetworkAnalysis/flashweave.jl


# subset mapping statistics to taxa of interest
keep = @pipe otudata|> combine(groupby(_,:tax_path), nrow=>:N_occur) |> filter(:N_occur=>o->o>20, _);
keep2 = @pipe otudata|> filter(:Middle_depth=>d->!(d in [1179.2,1283.3,966.1,900.5,350.5,10.5,0])) |> 
                        combine(groupby(_,:tax_path), :read_ani_mean=>minimum=>:min_ani) |> 
                        filter(:min_ani=>m->m>95, _);

tab = @pipe otudata |> filter(:tax_path=>p->p in keep.tax_path && p in keep2.tax_path, _);
unique(tab.tax_name)
filter!(:Middle_depth=>d->d!=900.5,tab);
filter!(:Middle_depth=>d->d<=1100.5,tab);



# prepare data for model
tab.ka = -1/1000 * tab.yrBP;
tab.bsANI = fill(0.0,nrow(tab));
tab2=filter(:Middle_depth=>d->34.5<=d<=300.5,tab);

for pth in unique(tab.tax_path)
    tmp = @pipe filter(:tax_path=>p->p==pth, tab) |> sort(_,:Middle_depth)
    tab[tab.tax_path.==pth,:bsANI] .= tmp[!,:read_ani_mean][end]
end
tab.dANI = tab.read_ani_mean .- tab.bsANI;

for pth in unique(tab2.tax_path)
    tmp = @pipe filter(:tax_path=>p->p==pth, tab2) |> sort(_,:Middle_depth)
    tab2[tab2.tax_path.==pth,:bsANI] .= tmp[!,:read_ani_mean][end]
end
tab2.dANI = tab2.read_ani_mean .- tab2.bsANI;



# mixed-effect model
using MixedModels, StatsModels
mod = fit(MixedModel, @formula(dANI ~ 1 + ka + + (1 + ka |tax_name)), filter(:Middle_depth=>d->d>152, tab))
mod = fit(MixedModel, @formula(dANI ~ 1 + ka + ka&ka + (1 + ka + ka&ka|tax_name)), filter(:Middle_depth=>d->d>=34.5, tab))

data=@pipe filter(:Middle_depth=>d->d>=34.5, tab)|>select(_,[:tax_name,:Middle_depth,:ka])
data.fitted=fitted(mod);
data.residual=residuals(mod);

# random effects
re = raneftables(mod)[1] |> DataFrame

# plot fit of random effect for each sp
for (i,sp) in enumerate(re.tax_name)
    tmp = @pipe filter(:tax_name=>n->n==sp,tab) |> sort(_,:Middle_depth)
    tb = filter(:tax_name=>n->n==sp,re)
    fig=Figure()
    Axis(fig[1,1], limits=(nothing,nothing,-3,3), title=sp)
    CairoMakie.lines!(Vector{Float64}(tmp.ka), tmp.dANI, linewidth=1.5, color=(:black,0.15))
    X = minimum(unique(tmp.ka)):0.01:maximum(unique(tmp.ka));
    Y = re[i,2].+re[i,3].*X .+ re[i,4].*X.*X;
    lines!(X,Y, linewidth=3,linestyle=:dot,color=:deepskyblue)
    save("./MicrobeProfiling/output/mixed_model_randef_$sp.pdf", fig);
end


# plot dmg for modelled sp
using CairoMakie
depths = tab.Middle_depth |> sort |> unique
fig=Figure(resolution=(420,1000))
Axis(fig[1, 1], title = "",limits=(0,0.4,-depths[end]-10,-depths[1]+150), yticks=-depths,
          xlabel="DNA damage", ylabel="depth (cm)", yreversed=false)
for dep in depths
    dt = @pipe filter(:Middle_depth => d->d==dep, tab); if nrow(dt)==0 continue end
    CairoMakie.hist!(Float64.(dt.damage); bins=0:0.005:0.4, normalization=:none, offset=-dep,
                     weights=fill(5,size(dt)[1]), color = (:teal, 0.5))
end
save("./MicrobeProfiling/output/hist_dmg_mixmodel_sp.pdf", fig);


# boxplot
fig=Figure(size=(650,420))
tk1=Float64.(sort(unique(tab.yrBP))[Not([2,3,4,6,8,10])]);
axs1=Axis(fig[1,1],limits=(nothing,nothing,-3,3), xlabel="Year BP", ylabel="ANI difference to 11.9 ka BP",
          xticks=(-tk1/1000, string.(round.(Int, tk1))), xticklabelrotation=π/2,
         xticklabelsize=12, yticklabelsize=12, xgridvisible=false, ygridvisible=false)

CairoMakie.boxplot!(axs1, tab.ka, tab.dANI, color=(:teal,0.6), width=0.15, strokewidth=0.5, strokecolor =:teal, show_outliers=false)
X = minimum(unique(tab.ka)):0.01:maximum(unique(tab.ka));
Y = coef(mod)[1].+coef(mod)[2].*X .+ coef(mod)[3].*X.*X;
lines!(axs1,X,Y, linewidth=2,color=(:navyblue,0.6))
save("./MicrobeProfiling/output/mixed_model_lp_ANI95_boxplot.pdf", fig);


# zoomed in boxplot
fig=Figure(size=(400,220))
tk1=Float64.(sort(unique(tab.yrBP)));
filtered_tab = @pipe filter(:Middle_depth=>n->n<150, tab)
axs1=Axis(fig[1,1],limits=(nothing,nothing,-3,3), xlabel="Year BP", ylabel="ANI difference to 11.9 ka BP",
          xticks=(-tk1/1000, string.(round.(Int, tk1))), xticklabelrotation=π/2,
         xticklabelsize=12, yticklabelsize=12, xgridvisible=false, ygridvisible=false)

CairoMakie.boxplot!(axs1, filtered_tab.ka, filtered_tab.dANI, color=(:teal,0.6), width=0.02, strokewidth=0.5, strokecolor =:teal, show_outliers=false)
X = minimum(unique(filtered_tab.ka)):0.01:maximum(unique(filtered_tab.ka));
Y = coef(mod)[1].+coef(mod)[2].*X .+ coef(mod)[3].*X.*X;
save("./MicrobeProfiling/output/mixed_model_lp_ANI95_boxplot_zoom.pdf", fig);



