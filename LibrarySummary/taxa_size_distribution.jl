#!/usr/bin/env julia
using Plots
using StatsPlots
using FileIO
using Measures
using DataFrames, CSV, JLD2, StatsPlots, Plots, Statistics

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw";add="_ANI92";
alg=".lca";
tag="$pip$alg$add"
frank="species" #genus, family

# Load the main DataFrame (containing: Label, n_reads, read_length_mean)
df = load_object("./InitialExploration/data/$tag.$frank.samplecontrol.jld2")

dropmissing!(df, [:n_reads, :read_length_mean])

# read taxa marked as contaminants in decontam (R)
contam1 = readlines("./deContamination/output/contam_taxa_name.txt")

# taxa in controls
ctrl= @pipe filter(:Middle_depth=>d->ismissing(d), df)
contam2 = ctrl.tax_name


# Load metadata
metadata = CSV.read("./metadata/HIBO_library_metadata.csv", DataFrame)

# Extract XX from Label and flag controls
metadata[!, :XX] = replace.(string.(metadata.Label), r"HB_" => "")
metadata[!, :is_control] = ismissing.(metadata.yrBP) .|| string.(metadata.yrBP) .== "NA"

# Create a mapping from "XX" to (Label, is_control)
label_map = Dict{String, Tuple{String, Bool}}()
for row in eachrow(metadata)
    xx = row.XX
    is_ctrl = row.is_control
    label = is_ctrl ? "$(row.Label) (control)" : row.Label
    label_map[xx] = (label, is_ctrl)
end

# Prepare list of plots
plots_list = Plots.Plot[]

# Loop through metadata in order
for row in eachrow(metadata)
    xx = row.XX
    label, is_control = label_map[xx]
    color = is_control ? :darkred : :darkblue

    # Filter df for this label
    df_sub = @pipe df[df.Label .== row.Label, :] 
    
    if nrow(df_sub) == 0
        @warn "All rows for $(row.Label) have missing values or were filtered out. Skipping plot."
        continue
    end
    # Create a density plot of n_reads, weighted by read_length_mean
    x = df_sub.read_length_mean
    weights = df_sub.n_reads

    p = density(x;
        weights=weights,
        trim=false,
        linewidth=0,
        fill=(0, :lightgray),
        linecolor=nothing,
        bandwidth = 0.5,
        title=label,
        grid=false,
        legend=false,
        xlabel="",
        ylabel="",
        xtickfont=font(8),   # Set x-axis tick font size
        ytickfont=font(8)    # Set y-axis tick font size

    )

    if !is_control
        df_sub_decon = @pipe df_sub |> filter(:tax_name=>n-> !(n in union(contam1,contam2)), _)
        x2 = df_sub_decon.read_length_mean
        w2 = df_sub_decon.n_reads

        df_sub_con = @pipe df_sub |> filter(:tax_name=>n-> n in union(contam1,contam2), _)
        x3 = df_sub_con.read_length_mean
        w3 = df_sub_con.n_reads

        density!(p, x2; weights=w2, bandwidth = 0.5, linewidth=1.5, linecolor=color)
        density!(p, x3; weights=w3, bandwidth = 0.5, linewidth=1.5, linecolor=color, linestyle=:dot)
    end

    push!(plots_list, p)
end

# Arrange plots in 11 rows x 4 columns (adjust as needed)
final_plot = plot(plots_list..., layout=(11, 4), left_margin=5mm, size=(1000, 1600), dpi=300)

# Save the final figure
savefig(final_plot, "./LibrarySummary/output/taxa_read_length_weighted.pdf")

