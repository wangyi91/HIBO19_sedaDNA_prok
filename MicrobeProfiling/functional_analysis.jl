#!/usr/bin/env julia
using JLD2
using DataFrames
using Pipe
using CSV
using DataFramesMeta, Plots
#using PlotlyKaleido, PlotlyBase
using LinearAlgebra

libs = CSV.File("./metadata/HIBO_library_metadata.csv"; missingstring="NA") |> DataFrame |> dropmissing;
Libs=["HB_32","HB_02","HB_33","HB_08","HB_38"] # Roman times
Libs_up=["HB_43","HB_25","HB_26","HB_27","HB_28","HB_01","HB_29","HB_30","HB_31","HB_32","HB_02","HB_33","HB_08","HB_38","HB_03"] # 407cm and above
Libs_low=setdiff(libs.Label,Libs_up)# below 407cm



# cat samples of interest, functional annotation only
data_combined = DataFrame();
for lib in Libs_up
    filepath = "../MAG/diamond/output/$(lib).anno.path.txt"

    # Check if file exists
    if !isfile(filepath)
        @warn "File not found, skipping: $filepath"
        continue
    end

    anno_path = CSV.File(filepath, header=0, delim='\t') |> DataFrame;
    apath = @pipe anno_path |> rename(_, [:contigname, :annopath]) |>
                DataFrames.transform(_, :annopath => ByRow(c -> split.(c, ';'; limit=3)) => [:annodb,:a1,:a2]);
    insertcols!(apath, :Label=>lib, 
                :Middle_depth=>libs[occursin.(lib, libs.Label), :Middle_depth][1], 
                :yrBP=>libs[occursin.(lib, libs.Label),:yrBP][1]) 
    data_combined = vcat(data_combined, apath, cols=:union)
end

data_comb = @pipe filter(:annodb=>d->d=="KEGG",data_combined) |> 
                  filter(:a1=>a->a=="Metabolism",_) |> 
                  DataFrames.transform(_, :a2=>ByRow(p->split(p,';'))=>[:a2,:a3,:a4,:i])

dt = @pipe data_comb |>
            combine(groupby(_, [:Label, :Middle_depth, :yrBP, :a3]), nrow=>"entries") |> 
            sort(_, [:a3, :Middle_depth]) |>
            unstack(_,:a3, :Label,:entries, fill=0);

# Compute row sums (assuming all columns are numeric)
row_sums = sum.(eachrow(dt[:, 2:end]))
dt.sum = row_sums

dt_sorted=@pipe sort(dt, :sum, rev=true)|> select(_, Not(:sum))

CSV.write("./MicrobeProfiling/output/table_gene_functions.csv",dt)

normby="sample" # normby="gene"
if normby=="sample"
    for col in names(dt_sorted)[2:end]
        dt_sorted[!, col] = (dt_sorted[!, col] .- mean(dt_sorted[!, col])) ./ std(dt_sorted[!, col])
    end
else # norm by gene, to plot
    cols = names(dt_sorted)[2:end]
    for col in cols
        dt_sorted[!, col] = Float64.(dt_sorted[!, col])
    end

    for i in 1:nrow(dt_sorted)
        row = collect(dt_sorted[i, cols])  # values to normalize
        norm_val = norm(row)
        if norm_val != 0
            dt_sorted[i, cols] .= row ./ norm_val
        end
    end
end

dt[:, :sum3] = dt.HB_02 .+ dt.HB_33 .+ dt.HB_08

if normby== "gene"
    gr(margins = 1.5Plots.cm, size=(600, 400))
    plot(heatmap(dt_sorted.a3, 
                 reverse(names(dt_sorted)[2:end]), 
                 reverse(dt_sorted[:, 2:end] |> Matrix |> transpose, dims=1)),
         xticks=nothing, size=(600, 400))#, fc=cgrad([:white,:dodgerblue4])))
    title!("normalised by $normby");
    savefig("./MicrobeProfiling/output/heatmap_function_407above_normby$normby.pdf")
end


### check what's special about HB_31 ###
df=DataFrame(mt_norm, names(dt_sorted[:,2:end-1]))
df.a3=dt_sorted.a3
sort!(df,:HB_31,rev=true)
### conclusion: these are from plants and fungus ###


