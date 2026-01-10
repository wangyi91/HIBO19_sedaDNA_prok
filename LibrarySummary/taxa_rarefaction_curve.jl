#!/usr/bin/env julia
using DataFrames, CSV, JLD2, StatsPlots, Plots, Statistics
using Random
using StatsFuns

include("./MicrobeProfiling/_get_annotated_df.jl")

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

using DataFrames
using StatsPlots
using Random

# Assume df is already loaded and has columns: Label, N_reads, tax_name

# Function to calculate rarefaction
function rarefaction_curve(counts::Vector{Int})
    n_total = sum(counts)
    curve = Float64[]
    for n in 1:1000:n_total
        # simulate sampling without replacement
        sampled = sample(counts, n; replace=false)
        push!(curve, length(unique(sampled)))
    end
    return curve
end

# Alternative, more efficient rarefaction using cumulative sum
function rarefaction_curve_efficient(counts::Vector{Int})
    # expand counts into a vector of tax_names repeated by counts
    taxa_expanded = repeat(1:length(counts), counts)
    # shuffle
    shuffle!(taxa_expanded)
    unique_counts = Int[]
    seen = Set{Int}()
    for x in taxa_expanded
        push!(seen, x)
        push!(unique_counts, length(seen))
    end
    return unique_counts
end

# Prepare data for plotting
sample_curves = Dict()
for s in unique(df.Label)
    subdf = df[df.Label .== s, :]
    counts = Vector{Int64}(subdf.N_reads)
    sample_curves[sample] = rarefaction_curve(counts)
end

# Plot all rarefaction curves
plt = plot(title="Rarefaction Curves", xlabel="Number of reads", ylabel="Observed taxa", legend=:topleft)
for (sample, curve) in sample_curves
    plot!(plt, 1:length(curve), curve, label=sample)
end

display(plt)


savefig("./LibrarySummary/output/rarefaction_curves.png")

