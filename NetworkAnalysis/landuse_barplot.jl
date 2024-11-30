#!/usr/bin/env julia

# This script make a simple barplot of land use data extracted from Ellis et al. and cited sources 
using DataFrames, JLD2
using Pipe
using CSV, CSVFiles
using CairoMakie

landuse = @pipe load("./metadata/land_use_anthromes.csv") |> DataFrame |> sort(_, :yrBP)

filter!(:yrBP=>y->!(y<0 && mod(y,10)!=0),landuse) # filter out modern data that are too dense

landuse.win = [0 for _=1:nrow(landuse)];
for i in 1:nrow(landuse)-1
    landuse.win[i] = landuse.yrBP[i+1]-landuse.yrBP[i]
end
landuse.win[nrow(landuse)] = landuse.win[nrow(landuse)-1];

long = stack(landuse, 2:ncol(landuse)-1);
long.stack = [findall(==(v), @pipe long.variable |> unique |> sort(_, rev=true))[1] for v in long.variable];
long.x = [findall(==(t), long.yrBP |> unique |> sort)[1] for t in long.yrBP];

xtl=unique(long[!,[:yrBP,:x]]); # a df for xticklabels

using DataStructures
landcolors= SortedDict("11"=>:firebrick,
                 "12"=>:firebrick1,
                 "23"=>:orchid4,
                 "24"=>:thistle2,
                 "32"=>:lightgoldenrod,
                 "33"=>:lightgoldenrodyellow,
                 "41"=>:orange,
                 "42"=>:goldenrod1,
                 "43"=>:wheat2,
                 "51"=>:forestgreen,
                 "52"=>:darkolivegreen3,
                 "53"=>:honeydew2,
                 "54"=>:wheat3,
                 "61"=>:lightblue1,
                 "62"=>:snow3)

landtypes = SortedDict("11"=>"Urban",
                 "12"=>"Mixed settlements",
                 "23"=>"Rainfed villages",
                 "24"=>"Pastoral villages",
                 "32"=>"Residential rainfed croplands",
                 "33"=>"Populated croplands",
                 "41"=>"Residential rangelands",
                 "42"=>"Populated rangelands",
                 "43"=>"Remote rangelands",
                 "51"=>"Residential woodlands",
                 "52"=>"Populated woodlands",
                 "53"=>"Remote woodlands",
                 "54"=>"Inhabited drylands",
                 "61"=>"Wild woodlands",
                 "62"=>"Wild drylands")

fig=Figure(resolution=(1000,500))
axs = Axis(fig[1, 1], limits=(minimum(xtl.x),maximum(xtl.x),0,1), xlabel="year BP", ylabel="%", xticks=unique(long.x), xreversed=true,
           xtickformat = tks->[ifelse(xtl[xtl.x.==tk,:].yrBP[1] in [-60,0,250,1950,11950], string(xtl[xtl.x.==tk,:].yrBP[1]),"") for tk in tks])

CairoMakie.barplot!(axs, long.x, long.value,
                    color=[landcolors[v] for v in long.variable], 
                    stack=long.stack, gap=0)

land_color = [PolyElement(color = clr, strokecolor = :transparent) for clr in collect(values(landcolors))];
Legend(fig[1,2], land_color, collect(values(landtypes)), "Landscape type")
save("./NetworkAnalysis/output/barplot_landuse.pdf", fig);




