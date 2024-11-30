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

# spaghetti plot or boxplot
fig=Figure(resolution=(740,420))
tk1=Float64.(sort(unique(tab.ka),rev=true)[Not([2,3,4,6,8,10])]);
axs1=Axis(fig[1,1],limits=(nothing,nothing,-3,3), title="Long-persisting species", xlabel="ka BP", ylabel="ANI difference to 11.9 ka BP",
         xticks=tk1, xtickformat=values-> [string(round(value*-1,digits=2)) for value in values], xticklabelrotation=π/2,
         xticklabelsize=12, yticklabelsize=12, xgridvisible=false, ygridvisible=false)

for (i,pth) in enumerate(unique(tab.tax_path))
    tmp = @pipe filter(:tax_path=>p->p==pth,tab) |> sort(_,:Middle_depth)
    #if !isempty(tmp[tmp.Middle_depth.==644.3,:dANI])&&tmp[tmp.Middle_depth.==300.5,:dANI][1]>1.5 clr=(:orangered,0.5) else clr=(:black, 0.15) end
    CairoMakie.lines!(axs1, Vector{Float64}(tmp.ka), tmp.dANI,
                      linewidth=1.5, color=(:black,0.15))
end
#boxplot!(axs1,tab.ka,tab.dANI,width=0.15,show_outliers=false)
#boxplot!(axs1,data.ka,data.residual,width=0.15,show_outliers=false)
X = minimum(unique(tab.ka)):0.01:maximum(unique(tab.ka));
Y = coef(mod)[1].+coef(mod)[2].*X .+ coef(mod)[3].*X.*X;
lines!(axs1,X,Y, linewidth=3,linestyle=:dot,color=:deepskyblue)
save("./MicrobeProfiling/output/mixed_model_lp_ANI95_spaghetti.pdf", fig);






















# some scripts for prelim exploration

using GLM, StatsBase
sts = DataFrame(tax_path=unique(tab.tax_path),
                   med_ani=fill(0.0,length(unique(tab.tax_path))),
                   min_ani=fill(0.0,length(unique(tab.tax_path))),
                   std_ani=fill(0.0,length(unique(tab.tax_path))),
                   Dani=fill(0.0,length(unique(tab.tax_path))),
                   slp=fill(0.0,length(unique(tab.tax_path))),
                   r2=fill(0.0,length(unique(tab.tax_path))))

for (i,pth) in enumerate(sts.tax_path)
    tmp = @pipe filter(:Middle_depth=>d->d>10.5,tab) |> filter(:tax_path=>p->p==pth,_)
    sts.med_ani[i] = median(tmp.read_ani_mean)
    sts.min_ani[i] = minimum(tmp.read_ani_mean)
    sts.Dani[i] = maximum(tmp.read_ani_mean)-minimum(tmp.read_ani_mean)
    sts.std_ani[i] = std(tmp.read_ani_mean)
    ols=lm(@formula(read_ani_mean~yrBP/1000),tmp)
    sts.slp[i]=coef(ols)[2]
    sts.r2[i]=r2(ols)
end
leftjoin!(tab, sts, on=:tax_path);


# plot stats of params for filtering
fig=Figure(resolution=(500,300))
axs1=Axis(fig[1,1],title="slp");axs2=Axis(fig[1,2],title="r2");
CairoMakie.hist!(axs1,unique(out.slp))
CairoMakie.hist!(axs2,unique(out.r2))
save("./MicrobeProfiling/output/hist_ani.pdf", fig);

# after examining the plots, further filter
data = @pipe tab |> filter(:r2=>r->r>0.45, _)
out=filter(:tax_name=>n->!(n in data.tax_name), tab); unique(out.tax_name)
out2 = @pipe out |> filter(:tax_name=>n->n in ["SIAJ01 sp009691845","SYLI01 sp009695465"],_)#,"SHWZ01 sp009693575"],_),"Methylomirabilis limnetica"],_)



using MixedModels, StatsModels
mod = fit(MixedModel, @formula(dANI ~ 1 + ka + ka&ka + (1 + ka + ka&ka|tax_name)), filter(:Middle_depth=>d->d>=34.5,data))


# set up plot
#=xtks=Float64.(sort(unique(tab.yrBP))[Not([2,3,4,6,8,10])])/1000;
fig = Figure(resolution=(500,700))
axs = Axis(fig[1,1], title="ΔANI(%) vs time", limits=(minimum(xtks)-0.25,nothing,nothing,nothing),
           xticks=xtks, xtickformat="{:.2f}", xticklabelsize=10, xticklabelrotation=π/2, yticks=-3:1:3,
           xgridwidth=0.5, ygridwidth=0.2)
ylabel = Label(fig[:, 0], ylab, rotation = pi/2, tellheight=false)
xlabel = Label(fig[2,:], "ka BP")=#

fig=Figure(resolution=(420,300))
axs=Axis(fig[1,1],limits=(nothing,nothing,-3,3))

for (i,pth) in enumerate(unique(data.tax_path))
    tmp = @pipe filter(:tax_path=>p->p==pth,data) |> sort(_,:Middle_depth)
    CairoMakie.lines!(axs, Vector{Float64}(tmp.ka), tmp.dANI,
                      linewidth=2, color=(:black, 0.15))
end
X=minimum(unique(data.ka)):0.01:maximum(unique(data.ka));
Y=coef(mod)[1].+coef(mod)[2].*X .+ coef(mod)[3].*X.*X;
lines!(X,Y,linewidth=2,linestyle=:dot,color=:royalblue)
save("./MicrobeProfiling/output/mixed_model.pdf", fig);


#outliers
fig=Figure(resolution=(420,300))
axs=[Axis(fig[1,i],limits=(nothing,nothing,-3,3)) for i in [1,2]]

for (i,pth) in enumerate(unique(out.tax_path))
    tmp = @pipe filter(:tax_path=>p->p==pth, out) |> sort(_,:Middle_depth)
    CairoMakie.lines!(axs[1], Vector{Float64}(tmp.ka), tmp.dANI,
                      linewidth=1.5, color=(:black, 0.15))
end
for (i,pth) in enumerate(unique(out2.tax_path))
    tmp = @pipe filter(:tax_path=>p->p==pth, out2) |> sort(_,:Middle_depth)
    CairoMakie.lines!(axs[2], Vector{Float64}(tmp.ka), tmp.dANI,
                      linewidth=1.5, color=(:black, 0.15))
end

save("./MicrobeProfiling/output/mixed_model_outliers.pdf", fig);












