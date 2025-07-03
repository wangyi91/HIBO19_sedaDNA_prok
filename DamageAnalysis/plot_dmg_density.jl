#!/usr/bin/env julia
using DataFrames
using JLD2
using CSV #CSVFiles
using Pipe
using Formatting

# plot damage as one-sided violin plot for all samples on a single figure

pip="amaw"; alg=".lca";add="_ANI92";frank="species";
pip="holi"; alg=".lca";add=".genus";frank="genus";

tag=pip*alg*add

df = load_object("./deContamination/data/$tag.$frank.jld2")

# choose dataframe to plot
depths = df.Middle_depth |> sort |> unique;
ages = df.yrBP |> sort |> unique;

colordic = Dict("Bacteria"=>:slateblue,"Archaea"=>:darkred, ""=>:teal,
                "Embryophyta"=>:olivedrab,"Chordata"=>:hotpink,"Arthropoda"=>:darkslategrey,
                "Non-plant eukaryote"=>:hotpink);

ftax="Archaea"; Nmin=50;
ftax="Bacteria"; Nmin=500;
#ftax="Embryophyta"
#ftax="Non-plant eukaryote"

using CairoMakie
# for holi data
fig = Figure(size=(500,1000))
Axis(fig[1, 1], limits=(0,0.45,-depths[end]-20,-depths[1]+70), yticks=-depths, 
     #title = "$ftax (after filtration)",
     title = "$ftax (reads / reference $frank > $Nmin)",
     xlabel="DNA damage", ylabel="depth (cm)", yreversed=false, yticklabelsize=12, xticklabelsize=12)
for dep in depths
    dt = @pipe filter(:Middle_depth => d->d==dep, df) |> filter(:tax_path=>p->occursin(ftax,p)) |> filter(:N_reads=>n->n>Nmin,_); 
    #dt = @pipe filter(:Middle_depth => d->d==dep, df) |> filter(:tax_path=>p->occursin("Eukaryota",p)) |> filter(:tax_path=>p->!occursin("Embryophyta",p))
    if nrow(dt)==0 continue end
    #w = Float64.(sqrt.(dt.N_reads))/10
    w = Float64.(dt.N_reads)
    #CairoMakie.hist!(Float64.(dt.damage); bins=0:0.005:0.4, normalization=:none, offset=-dep, #weights=w, #fill(3,size(dt)[1]), 
                     #color = (:firebrick, 0.5))
    CairoMakie.density!(Float64.(dt.damage), npoints = 5000, offset=-dep, weights = w,
                        color = (colordic[ftax], 0.4), bandwidth = 0.002)
end
save("./DmgMixtureModel/output/dmg_density_$(tag)_$(ftax)_afterfilt.pdf", fig);

# CairoMakie.scatter!(stats.lower, stats.Middle_depth, color = (:red, 0.4), markersize = 8, marker = '▲')
# CairoMakie.lines!(0..0.3, x -> 1/GLM.coef(reg)[2]*x - GLM.coef(reg)[1]/GLM.coef(reg)[2], color = :red, linewidth = 2) # reg is done with yrBP not depth. better add line manually. 

# plot by age not depth
fig = Figure(size=(500,1000))
Axis(fig[1, 1], limits=(0,0.45,-ages[end]-200,-ages[1]+700), yticks=-ages, ytickformat=tks->[string(round(Int,tk)) for tk in tks],
     #title = "$ftax (after filtration)",
     title = "$ftax (reads per reference $frank > $Nmin)",
     xlabel="DNA damage", ylabel="year BP", yreversed=false, yticklabelsize=12, xticklabelsize=12)
for age in ages
    dt = @pipe filter(:yrBP => a->a==age, df) |> filter(:tax_path=>p->occursin(ftax,p)) |> filter(:N_reads=>n->n>Nmin,_);
    if nrow(dt)==0 continue end
    #w = Float64.(dt.N_reads)/200# weights for bacteria
    w = Float64.(dt.N_reads)/30# weights for archaea
    #w = Float64.(sqrt.(dt.N_reads))/2 # weights for plants
    CairoMakie.hist!(Float64.(dt.damage); bins=0:0.004:0.4, normalization=:none, offset=-age, weights=w, color = (colordic[ftax], 0.4))#bins=0:0.004:0.4 for bacteria
    #CairoMakie.density!(Float64.(dt.damage), npoints = 5000, offset=-age, weights = w, color = (colordic[ftax], 0.4), bandwidth = 5e-4)
end
save("./DamageAnalysis/output/dmg_hist_$(tag)_$(ftax)_age.pdf", fig);



# plot a partial plot by age not depth
ruler=[0,300,600,900,1200,1500,1800,2100,2400,2700,3000,3300,3600,3900]

age_dict = Dict(ages[i] => ruler[i] for i in 1:14)

fig = Figure(size=(250,600))
Axis(fig[1, 1], limits=(0,0.2,-ruler[end]-50,-ruler[1]+300), yticks=-ruler, ytickformat=tks->[string(round(Int,tk)) for tk in tks],
     #title = "$ftax (after filtration)",
     title = "$ftax (reads per reference $frank > $Nmin)",
     xlabel="DNA damage", ylabel="year BP", yreversed=false, yticklabelsize=12, xticklabelsize=12)
for age in ages[1:14]
    dt = @pipe filter(:yrBP => a->a==age, df) |> filter(:tax_path=>p->occursin(ftax,p)) |> filter(:N_reads=>n->n>Nmin,_);
    if nrow(dt)==0 continue end
    #w = Float64.(dt.N_reads)/400# weights for bacteria
    w = Float64.(dt.N_reads)/30# weights for archaea
    #w = Float64.(sqrt.(dt.N_reads))/2 # weights for plants
    CairoMakie.hist!(Float64.(dt.damage); bins=0:0.004:0.2, normalization=:none, offset=-age_dict[age], weights=w, color = (colordic[ftax], 0.6))#bins=0:0.004:0.4 for bacteria
    #CairoMakie.density!(Float64.(dt.damage), npoints = 5000, offset=-age, weights = w, color = (colordic[ftax], 0.4))
end
save("./DamageAnalysis/output/dmg_hist_$(tag)_$(ftax)_age_partial.pdf", fig);

