#!/usr/bin/env julia
using DataFrames, JLD2, Pipe, CSV, StatsBase

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw";add="_ANI92";
alg=".lca";
tag="$pip$alg$add"
frank="species" #genus, family


# load datatable that contain both samples and controls
df=load_object("./InitialExploration/data/$tag.$frank.samplecontrol.jld2")

# read taxa marked as contaminants in decontam (R)
contam1 = readlines("./deContamination/output/contam_taxa_name.txt")

# taxa in controls
ctrl= @pipe filter(:Middle_depth=>d->ismissing(d), df)
contam2 = ctrl.tax_name

# remove taxa in controls and in "contam" list, drop controls
df_clean = @pipe df |> filter(:tax_name=>n->!(n in union(contam1,contam2)), _) |> 
                       filter(:Label => l -> !(l in ctrl.Label),_)


# compare samples before and after decontamination, in PCA plot
using Statistics
using MultivariateStats
using CairoMakie
using CategoricalArrays

# 1. Pivot to sample x tax matrix for each dataframe 
function pivot_tax_table(df::DataFrame)
    unstack(df, :Label, :tax_name, :N_reads, fill =0)
end

df_pivot = @pipe pivot_tax_table(df) |> filter(:Label => l -> l in df_clean.Label, _)
df_clean_pivot = pivot_tax_table(df_clean)

# 2. Align tax columns (fill missing taxa with 0s) 
for col in setdiff(names(df_pivot), names(df_clean_pivot))
    df_clean_pivot[!, col] = fill(0, nrow(df_clean_pivot))
end



# 3. Extract sample names and prepare data matrix 
group_shape = vcat(fill(:cross, length(df_pivot.Label)), fill(:utriangle, length(df_pivot.Label)))
group_color = vcat(unique(df_clean.yrBP),unique(df_clean.yrBP))

X = vcat(df_pivot[:, Not(:Label)], df_clean_pivot[:, Not(:Label)]) |> Matrix |> transpose

# 4. Optional: log transform or normalize 
X_log = log1p.(X')



# 5. Perform PCA 
M = fit(PCA, X_log; maxoutdim=2)
pj = @pipe projection(M) |> DataFrame(_, :auto);

X_pca = MultivariateStats.transform(M, X_log)


# 6. Plot 
using Plots
gr(margins = 0.4Plots.cm, size=(400, 325));

p=Plots.plot(X_pca[1,:], X_pca[2,:], seriestype = :scatter,
             group = vcat(df_pivot.Label,df_pivot.Label),
             markershape = group_shape,
             marker_z=group_color,
             markersize = 4,              
             alpha = 0.9,                
             markerstrokewidth = 0,     
             grid = false,             
             xlabel = "PC1",          
             ylabel = "PC2",         
             label="", legend=true);
savefig("./deContamination/output/pca_decontam.pdf");


# 7. Save filtered dataset
save_object("./deContamination/data/$tag.$frank.jld2", df_clean)









