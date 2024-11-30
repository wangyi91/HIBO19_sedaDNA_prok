#!/usr/bin/env julia
using Interpolations
using CSV, CSVFiles
using DataFrames
using Pipe: @pipe

temp12k = @pipe load("./metadata/temp12k_allmethods_percentiles.csv") |> DataFrame |> select(_, [:ages, Symbol("30N_to_60N_median")]);
temp1850= @pipe DataFrame(CSV.File("./metadata/NCEI_NH_Land_AverageTemp.csv", delim = ',', skipto=6, header=5))|>
    DataFrames.transform(_, :Date=>ByRow(d->round((195000-d)/100)-(mod(d,100)-1)/12)=>:yrBP)|>
    DataFrames.transform(_, :Anomaly=>ByRow(n->n+temp12k[1,2])=>:Anomaly_adj) # change to relative to 1800-1900 median (before is 1901-2000 average)

xrf = @pipe load("./metadata/2021-040_Harms-et-al_3_XRF_Element_Data_Composite_Core.csv") |> DataFrame
sic = @pipe load("./metadata/2021-040_Harms-et-al_4_Silica(bSi)_Carbon(TiC_ToC)_Combined_Section_data.csv") |> DataFrame |> 
            rename(_, Dict("bSi_Content_[wt%]           " =>"bSi_Content_[wt%]"))
depcorr= @pipe load("./metadata/depth_event_correction.csv") |> DataFrame
agedepthmodel = @pipe load("./metadata/age_depth_model_specification.csv") |> DataFrame |> filter(:yrBP=>y->!ismissing(y),_)|> sort(_,:yrBP) |> select(_, [:z,:yrBP,:sigma])

# add corrected depth for env variables
xrf_corr = @pipe xrf |>
    DataFrames.transform(_, "Comp Depth [cm]" => ByRow(d->linear_interpolation(depcorr.Comp_depth_cm, depcorr.Comp_depth_corrected, extrapolation_bc=Line())(d)) => :depth_corrected)|> select(_,[:depth_corrected,:S])

sic_corr = @pipe sic |>
    DataFrames.transform(_, "Comp_Depth_[cm]" => ByRow(d->linear_interpolation(depcorr.Comp_depth_cm, depcorr.Comp_depth_corrected, extrapolation_bc=Line())(d)) => :depth_corrected)|> select(_,["depth_corrected","bSi_Content_[wt%]","TOC_[wt%]","TIC_[wt%]"])

temp12k_corr = @pipe temp12k |>
    DataFrames.transform(_, :ages=>ByRow(a->linear_interpolation(agedepthmodel.yrBP, agedepthmodel.z, extrapolation_bc=Line())(a)) => :depth_corrected)
temp1850_corr = @pipe temp1850 |>
    DataFrames.transform(_, :yrBP=>ByRow(a->linear_interpolation(agedepthmodel.yrBP, agedepthmodel.z, extrapolation_bc=Line())(a)) => :depth_corrected)


using CairoMakie
xmin=[0,0,0,nothing]
xmax=[1000,nothing,nothing,nothing]
fig=Figure()
resize!(fig, 450,900)
axs=[Axis(fig[1,i],limits=(xmin[i],xmax[i],-1330,0), xgridvisible=false, leftspinevisible=false, rightspinevisible=false) for i in 1:4]
hideydecorations!.(axs)
colgap!(fig.layout, 10)
lines!(axs[1],xrf_corr.S,xrf_corr.depth_corrected*(-1),color=:orange, linewidth=0.8)
lines!(axs[2],sic_corr."bSi_Content_[wt%]",sic_corr.depth_corrected*(-1), color=:dodgerblue, linewidth=1)
lines!(axs[3],sic_corr."TOC_[wt%]",sic_corr.depth_corrected*(-1), color=:black, linewidth=1)
lines!(axs[4],temp1850_corr.Anomaly_adj,temp1850_corr.depth_corrected*(-1), color=:salmon, linewidth=1)
lines!(axs[4],temp12k_corr."30N_to_60N_median",temp12k_corr.depth_corrected*(-1), color=:darkred, linewidth=1)

save("./NetworkAnalysis/output/envs.pdf",fig)
