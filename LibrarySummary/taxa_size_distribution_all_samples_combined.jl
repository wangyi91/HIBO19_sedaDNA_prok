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


# === PREPARE LEGEND PLOT AND COMBINED PLOT ===

# Combined data arrays for the main plot
x_all_noncontam = Float64[]
w_all_noncontam = Float64[]
x_all_contam = Float64[]
w_all_contam = Float64[]

# Loop again to collect data for combined plot
for row in eachrow(metadata)
    xx = row.XX
    label, is_control = label_map[xx]
    color = is_control ? :darkred : :darkblue

    df_sub = df[df.Label .== row.Label, :]
    if nrow(df_sub) == 0
        continue
    end

    df_sub_decon = filter(:tax_name => n -> !(n in union(contam1, contam2)), df_sub)
    df_sub_con = filter(:tax_name => n -> n in union(contam1, contam2), df_sub)

    append!(x_all_noncontam, df_sub_decon.read_length_mean)
    append!(w_all_noncontam, df_sub_decon.n_reads)

    append!(x_all_contam, df_sub_con.read_length_mean)
    append!(w_all_contam, df_sub_con.n_reads)
end

# === LEGEND PLOT ===
# This is a dummy plot just to hold the legend
p_legend = plot(
    legend=:outertopright,
    title="Legend",
    grid=false,
    xlabel="",
    ylabel="",
    xtickfont=font(8),
    ytickfont=font(8)
)

# Dummy data for legend
plot!(p_legend, [NaN], label="Sample (non-contaminant)", linecolor=:darkblue, linewidth=1.5)
plot!(p_legend, [NaN], label="Control (non-contaminant)", linecolor=:darkred, linewidth=1.5)
plot!(p_legend, [NaN], label="Sample (contaminant)", linecolor=:darkblue, linewidth=1.5, linestyle=:dot)
plot!(p_legend, [NaN], label="Control (contaminant)", linecolor=:darkred, linewidth=1.5, linestyle=:dot)

# === COMBINED PLOT FOR ALL SAMPLES ===
combined_plot = plot(
    title="All Samples Combined",
    xlabel="Read Length Mean",
    ylabel="Density",
    legend=:topright,
    grid=false,
    xtickfont=font(8),
    ytickfont=font(8),
    size=(800, 300)
)

# Plot both contaminant and non-contaminant distributions
density!(combined_plot, x_all_noncontam; weights=w_all_noncontam, linewidth=1.5, linecolor=:darkblue, label="Sample (non-contaminant)")
density!(combined_plot, x_all_contam; weights=w_all_contam, linewidth=1.5, linecolor=:darkblue, linestyle=:dot, label="Sample (contaminant)")

# Optional: split control from sample if you want 4 lines in legend
# For clarity, above we just aggregate them together

# === COMBINE ALL PLOTS TOGETHER ===

# Append the legend and main combined plot to your plots list
# For example, make a 12x4 layout: 11 rows for samples, row 12 for combined + legend

# First make a blank plot to occupy grid positions
blank = plot()

# New layout: 12 rows, 4 columns
plots_list_with_combined = copy(plots_list)
push!(plots_list_with_combined, combined_plot)
push!(plots_list_with_combined, p_legend)
push!(plots_list_with_combined, blank)
push!(plots_list_with_combined, blank)

final_plot = plot(plots_list_with_combined...,
    layout=(12, 4),
    left_margin=5mm,
    size=(1000, 1800),
    dpi=300
)

# Save
savefig(final_plot, "./LibrarySummary/output/taxa_read_length_weighted_combined.pdf")

