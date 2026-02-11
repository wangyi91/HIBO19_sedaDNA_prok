#!/usr/bin/env julia
using DataFrames
using JLD2
using CSV #CSVFiles
using Pipe
using Formatting

# plot damage as offset histogram for all samples on a single figure

pip="amaw"; alg=".lca";add="_ANI92";frank="species";
tag=pip*alg*add

df = load_object("./deContamination/data/$tag.$frank.jld2")

# choose dataframe to plot
depths = df.Middle_depth |> sort |> unique;
ages = df.yrBP |> sort |> unique;

colordic = Dict("Bacteria"=>:royalblue3,"Archaea"=>:deeppink4, ""=>:teal,
                "Embryophyta"=>:olivedrab,"Chordata"=>:hotpink,"Arthropoda"=>:darkslategrey,
                "Non-plant eukaryote"=>:hotpink);

ftax="Archaea"; Nmin=50;scale_reads = 10_000
ftax="Bacteria"; Nmin=500;scale_reads = 50_000

# weight for histogram
weightdict = Dict(
    "Archaea"     => x -> x / 30,
    "Bacteria"    => x -> x / 200,
    "Embryophyta" => x -> sqrt(x) / 2
)


using CairoMakie
# plot by age not depth
fig = Figure(size=(370,1000))
ax = Axis(fig[1, 1], limits=(0,0.35,-ages[end]-200,-ages[1]+700), yticks=(-ages, string.(round.(Int, ages))), 
     title = "$ftax (reads per reference $frank > $Nmin)", titlefont = :regular, titlesize = 14,
     xlabel="DNA damage", ylabel="Year BP", yreversed=false, yticklabelsize=12, xticklabelsize=12, xlabelsize=13, ylabelsize=13)
     
ax.xgridvisible = false
for age in ages
    dt = @pipe filter(:yrBP => a->a==age, df) |> filter(:tax_path=>p->occursin(ftax,p)) |> filter(:N_reads=>n->n>Nmin,_);
    if nrow(dt)==0 continue end
    w = weightdict[ftax].(Float64.(dt.N_reads))
    CairoMakie.hist!(Float64.(dt.damage); bins=0:0.004:0.4, normalization=:none, offset=-age, weights=w, color = (colordic[ftax], 0.4))
end

# --- scale bar ---
scale_height = weightdict[ftax](scale_reads)

# position (top-right corner, tweak if needed)
xpos = 0.25
y0   = -ages[1] + 200   # small offset from top

lines!(ax,
       [xpos, xpos],
       [y0, y0 + scale_height],
       color = (colordic[ftax], 0.4),
       linewidth = 5)

text!(ax,
      xpos + 0.01,
      y0 + scale_height / 2,
      text = "$scale_reads reads",
      align = (:left, :center),
      fontsize = 11)

save("./DamageAnalysis/output/dmg_hist_$(tag)_$(ftax)_age.pdf", fig);



# --- partial plot for top samples only, evenly spaced ---
# weight for histogram
weightdict = Dict(
    "Archaea"     => x -> x / 45,
    "Bacteria"    => x -> x / 300
)

ruler=[0,300,600,900,1200,1500,1800,2100,2400,2700,3000,3300,3600,3900]

age_dict = Dict(ages[i] => ruler[i] for i in 1:14)

fig = Figure(size=(250,600))
ax=Axis(fig[1, 1], limits=(0,0.2,-ruler[end]-50,-ruler[1]+500), yticks=(-ruler, string.(round.(Int, ages[1:14]))),
     title = "$ftax, top 14 samples", titlesize=14, titlefont = :regular,
     xlabel="DNA damage", ylabel="Year BP", yreversed=false, yticklabelsize=12, xticklabelsize=12, ylabelsize=12, xlabelsize=12)
ax.xgridvisible = false
for age in ages[1:14]
    dt = @pipe filter(:yrBP => a->a==age, df) |> filter(:tax_path=>p->occursin(ftax,p)) |> filter(:N_reads=>n->n>Nmin,_);
    if nrow(dt)==0 continue end
    w = weightdict[ftax].(Float64.(dt.N_reads))
    CairoMakie.hist!(Float64.(dt.damage); bins=0:0.004:0.2, normalization=:none, offset=-age_dict[age], weights=w, color = (colordic[ftax], 0.4))
end

# --- scale bar ---
scale_height = weightdict[ftax](scale_reads)

# position (top-right corner, tweak if needed)
xpos = 0.1
y0   = -ages[1] + 50   # small offset from top

lines!(ax,
       [xpos, xpos],
       [y0, y0 + scale_height],
       color = (colordic[ftax], 0.4),
       linewidth = 5)

text!(ax,
      xpos + 0.01,
      y0 + scale_height / 2,
      text = "$scale_reads reads",
      align = (:left, :center),
      fontsize = 11)

save("./DamageAnalysis/output/dmg_hist_$(tag)_$(ftax)_age_partial.pdf", fig);


