#!/usr/bin/env julia
using Plots
default(fontfamily = "DejaVu Sans")
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

# Combined data arrays for the main plot
x_all = df.read_length_mean
w_all = df.n_reads

control_labels = metadata.Label[metadata.is_control .== true]
df_controls = filter(:Label => label -> label in control_labels, df)
x_ctrl = df_controls.read_length_mean
w_ctrl = df_controls.n_reads

df_cont = filter(:tax_name=>n-> n in union(contam1,contam2), df)
df_noncont = filter(:tax_name=>n-> !(n in union(contam1,contam2)), df)
x_noncontam = df_noncont.read_length_mean
w_noncontam = df_noncont.n_reads
x_contam = df_cont.read_length_mean
w_contam = df_cont.n_reads

combined_plot = plot(
    title="",
    xlabel="Average Read Length",
    ylabel="Density",
    legend=:topright,
    grid=false,
    xtickfont=font(8),
    ytickfont=font(8),
    left_margin=5mm, bottom_margin=5mm,
    size=(400, 300)
)

# Plot both contaminant and non-contaminant distributions
density!(combined_plot, x_all; weights=w_all, linewidth=0, fill=(0,:darkgrey),fillalpha=0.5, linecolor=nothing, bandwidth = 0.5, label="All taxa")
density!(combined_plot, x_ctrl; weights=w_ctrl, linewidth=0, fill=(0,:rosybrown),fillalpha=0.5, linecolor=nothing, bandwidth = 0.5, label="Taxa in controls")
density!(combined_plot, x_contam; weights=w_contam, linewidth=1.5, linecolor=:darkred, linestyle=:dot, bandwidth = 0.5, label="Identified contaminants")
density!(combined_plot, x_noncontam; weights=w_noncontam, linewidth=1.5, linecolor=:darkblue, bandwidth = 0.5, label="Non-contaminants")

# Save
savefig(combined_plot, "./LibrarySummary/output/taxa_read_length_weighted_combined.pdf")

