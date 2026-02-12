#!/usr/bin/env julia
using Interpolations
using CSV, CSVFiles
using DataFrames, FlashWeave
using Pipe: @pipe

function compile_metadata()
    temp12k = @pipe load("./metadata/temp12k_allmethods_percentiles.csv") |> DataFrame |> select(_, [:ages, Symbol("30N_to_60N_median")]);
    samples = @pipe load("./metadata/HIBO_library_metadata.csv") |> DataFrame |> dropmissing |> select(_, [:Label, :Middle_depth]) #|>
        #DataFrames.transform(_, :Top_depth=>ByRow(t->t+0.5)=>:Middle_depth_raw)
    dates = @pipe load("./metadata/HIBO_library_metadata_dated.csv") |> DataFrame |> dropmissing |> select(_, [:Label, :Middle_depth, :yrBP, :sigma])
    xrf = @pipe load("./metadata/2021-040_Harms-et-al_3_XRF_Element_Data_Composite_Core.csv") |> DataFrame
    sic = @pipe load("./metadata/2021-040_Harms-et-al_4_Silica(bSi)_Carbon(TiC_ToC)_Combined_Section_data.csv") |> DataFrame |>
        rename(_, Dict("bSi_Content_[wt%]           " =>"bSi_Content_[wt%]"))

    landuse = @pipe load("./metadata/land_use_anthromes.csv") |> DataFrame |> sort(_, :yrBP)
    libs = load("./metadata/HIBO_library_metadata_dated.csv") |> DataFrame |> dropmissing

    # add geochemical data by linear interpolation
    meta = @pipe samples |>
        DataFrames.transform(_, :Middle_depth => ByRow(d -> linear_interpolation(xrf."Comp Depth [cm]",xrf.S)(d)) => :S) |>
        DataFrames.transform(_, :Middle_depth => ByRow(d -> linear_interpolation(sic."Comp_Depth_[cm]",sic."bSi_Content_[wt%]",extrapolation_bc=Line())(d)) => :bSi) |>
        DataFrames.transform(_, :Middle_depth => ByRow(d -> linear_interpolation(sic."Comp_Depth_[cm]",sic."TOC_[wt%]",extrapolation_bc=Line())(d)) => :TOC) |>
        DataFrames.transform(_, :Middle_depth => ByRow(d -> linear_interpolation(sic."Comp_Depth_[cm]",sic."TIC_[wt%]",extrapolation_bc=Line())(d)) => :TIC)

    # join shifted depths (dated) and discard raw depths
    #meta = @pipe meta |> select(_, Not([:Top_depth,:Middle_depth_raw])) |> leftjoin(_, dates, on=:Label) |> leftjoin(_, select(libs, [:Label,:period]), on=:Label)
    meta = @pipe meta |>select(_, Not(:Middle_depth)) |> leftjoin(_, dates, on=:Label) |> leftjoin(_, select(libs, [:Label,:period]), on=:Label)

    # add land use data by interpolation on date
    for var in names(landuse)[2:end]
        meta = @pipe meta |> DataFrames.transform(_, :yrBP => ByRow(y -> linear_interpolation(landuse.yrBP,landuse[!,var],extrapolation_bc=Line())(y)) => Symbol(var))
    end

    # add temperature data by interpolation on date
    meta = @pipe meta |> DataFrames.transform(_, :yrBP => ByRow(y -> linear_interpolation(temp12k.ages,temp12k[!,"30N_to_60N_median"],extrapolation_bc=Line())(y)) => :N30_to_60_median); meta.N30_to_60_median[1:3]=[1,1.4,0] # correct temp after 0 year BP
    
    select!(meta, 1:8,:N30_to_60_median,:)
    return(meta)
end

