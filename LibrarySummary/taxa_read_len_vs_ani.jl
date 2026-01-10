#!/usr/bin/env julia
using Plots
default(fontfamily = "DejaVu Sans")
using StatsPlots
using FileIO
using Measures
using DataFrames, CSV, JLD2, StatsPlots, Plots, Statistics
using CairoMakie

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw";add="_ANI92";
alg=".lca";
tag="$pip$alg$add"
frank="species" #genus, family

# Load the main DataFrame (containing: Label, n_reads, read_length_mean)
df = load_object("./deContamination/data/$tag.$frank.jld2")
dropmissing!(df, [:n_reads, :read_length_mean])

# Load metadata
metadata = CSV.read("./metadata/HIBO_library_metadata.csv", DataFrame)

# Extract XX from Label and flag controls
metadata[!, :XX] = replace.(string.(metadata.Label), r"HB_" => "")
metadata[!, :is_control] = ismissing.(metadata.yrBP) .|| string.(metadata.yrBP) .== "NA"


# Combined data arrays for the main plot
x_all = df.damage
#x_all = df.read_length_mean
y_all = df.n_reads
z = df.read_ani_mean


z_ranges = [(92, 94), (94, 96), (96, 98), (98, 100)]

# Create a figure 
f = Figure(size = (1000, 800))

axes = []
for (i, (zmin, zmax)) in enumerate(z_ranges)
    row, col = Tuple(CartesianIndices((2, 2))[i])
    ax = Axis(f[row, col]; 
              yscale=log10, 
              #xlabel="Read Length Mean", 
              xlabel="DNA Damage", 
              ylabel="Number of Reads", 
              title = "ANI $zmin - $zmax",
             limits=(minimum(x_all), maximum(x_all), minimum(y_all), maximum(y_all)))
    mask = (z .>= zmin) .& (z .< zmax)
    CairoMakie.scatter!(ax, x_all[mask], y_all[mask];
             color = Float64.(z[mask]),   # explicitly Float64
             colormap = cgrad(:viridis, rev=true),
             colorrange = (92, 100),
             markersize = 6, alpha = 0.6)
    push!(axes, ax)
end

# Add a shared colorbar (take the first scatterplot's colormap)
Colorbar(f[:, end+1], axes[1].scene.plots[1];
         label = "ANI", vertical = true)

# Save
save("./LibrarySummary/output/taxa_dmg_vs_ani_panels.png", f)
