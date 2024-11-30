#!/usr/bin/env julia
using DataFrames, JLD2
using Pipe: @pipe

# This script is to model damage vs time and select for living species

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

alg=".lca";pip="amaw";add="_ANI92";tag=pip*alg*add
otudata = load_object("./InitialExploration/data/$tag.jld2")

keep = @pipe otudata|> combine(groupby(_,:tax_path), nrow=>:N_occur) |> filter(:N_occur=>o->o>=6, _);
dt = @pipe otudata |> filter(:tax_path=>p->p in keep.tax_path)
dt.yr=Float64.(dt.yrBP).+62.83
dt.dmg_adj=Float64.(dt.damage).-0.005


using GLM, StatsBase
fits = DataFrame(tax_path=unique(dt.tax_path),
                   slp=fill(0.0,length(unique(dt.tax_path))),
                   r2=fill(0.0,length(unique(dt.tax_path))));

for (i,pth) in enumerate(fits.tax_path)
    tmp = @pipe filter(:tax_path=>p->p==pth, dt)
    ols=lm(@formula(dmg_adj ~ 0 + yr),tmp)
    fits.slp[i]=coef(ols)[1]
    fits.r2[i]=r2(ols)
end

# visual inspection
using CairoMakie
fig=Figure(resolution=(800,400))
ax1=Axis(fig[1,1],xlabel="Damage accumulation rate (year⁻¹)",ylabel="R²")
ax2=Axis(fig[1,2],xlabel="Damage accumulation rate (year⁻¹)")
scatter!(ax1, fits.slp, fits.r2, markersize=3)
hist!(ax2, fits.slp, bins=100)
vlines!(ax1,0.65e-5, color=:darkred,linewidth=0.5, linestyle=:dash)
vlines!(ax2,0.65e-5, color=:darkred,linewidth=0.5, linestyle=:dash)
save("./DmgMixtureModel/output/dmg_lm_stats.pdf",fig)



# define and output damage of living taxa
# by linear increase of damage for holocene taxa
lv1 = @pipe filter(:slp=>s->s<0.65e-5,fits)|> filter(:r2=>r->r>0.4,_).tax_path

# by number of occurrences
lv2 = @pipe otudata |> combine(groupby(_,:tax_path), nrow=>:noccur)|>filter(:noccur=>n->n>8,_).tax_path

# Pleistocene taxa
lv3 = @pipe otudata|> combine(groupby(_,:tax_path), nrow=>:N_occur, :Middle_depth=>minimum=>:mindepth,:damage=>maximum=>:maxdmg) |>
            filter([:N_occur,:mindepth,:maxdmg]=>(o,t,d)->(o>=1&&t>1198&&d<0.05), _).tax_path

liv = filter(:tax_path=>p->p in [intersect(lv1, lv2);lv3], otudata) 

save_object("./DmgMixtureModel/output/liv.jld2",liv)
liv=load_object("./DmgMixtureModel/output/liv.jld2")
wide_dmg = @pipe liv |> unstack(_, :Middle_depth, :tax_name, :damage, fill=0)
CSV.write("./DmgMixtureModel/output/liv_dmg.csv", wide_dmg)

# send to local then use R.lda/
rsync -avP tvg137@dandycomp01fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/DmgMixtureModel/output/liv_*.csv /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.lda/Data/


# use R.lda/heatmap_dmg.R to plot and cut into 4 groups
# or, run here
mt = @pipe wide_dmg[:,2:end] |> (_.>0).*0.05.+_ |> Matrix # adjust dmg by adding 0.05 so that it's more separated from 0
dist = pairwise(euclidean, mt, dims=2)
ht = hclust(dist, linkage=:complete)

















# plot dmg profile of living species
fig = Figure(resolution=(420,950))
depths = otudata.Middle_depth |> sort |> unique
Axis(fig[1, 1], title="", limits = (0,0.4,-depths[end]-10,-depths[1]+50), 
     yticks=-depths, xlabel="DNA damage", ylabel="depth (cm)", yreversed=false)
for dep in depths
    dt = @pipe filter(:Middle_depth => d->d==dep, liv); if nrow(dt)==0 continue end
    CairoMakie.hist!(Float64.(dt.damage); bins=0:0.005:0.4, normalization=:none, offset=-dep,
                     weights=fill(5,size(dt)[1]), color = (:teal, 0.5))
end
save("./DmgMixtureModel/output/hist_dmg_liv.pdf", fig);













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
mod = fit(MixedModel, @formula(dANI ~ 1 + ka + ka&ka + (1 + ka + ka&ka|tax_name)), filter(:Middle_depth=>d->d>=34.5, tab))

data=@pipe filter(:Middle_depth=>d->d>=34.5, tab)|>select(_,[:tax_name,:Middle_depth,:ka])
data.fitted=fitted(mod);
data.residual=residuals(mod);


# plot
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
boxplot!(axs1,tab.ka,tab.dANI,width=0.15,show_outliers=false)
#boxplot!(axs1,data.ka,data.residual,width=0.15,show_outliers=false)
X = minimum(unique(tab.ka)):0.01:maximum(unique(tab.ka));
Y = coef(mod)[1].+coef(mod)[2].*X .+ coef(mod)[3].*X.*X;
lines!(axs1,X,Y, linewidth=3,linestyle=:dot,color=:deepskyblue)
save("./MicrobeProfiling/output/mixed_model_noccur20_resid.pdf", fig);






















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












