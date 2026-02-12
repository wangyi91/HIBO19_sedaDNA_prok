#!/usr/bin/env julia
using DataFrames, JLD2, Pipe, CSV, StatsBase

include("./MicrobeProfiling/_get_annotated_df.jl")
include("./NetworkAnalysis/_graph_process.jl")

pip="amaw";add="_ANI92";
alg=".lca";
tag="$pip$alg$add"
rank="species"


# Step 1: load datatable that contain both samples and controls
df = load_object("./InitialExploration/data/$tag.$rank.samplecontrol.jld2")


# Group by sample (Label) and taxon
grouped = combine(groupby(df, [:Label, :tax_name]), :N_reads => sum => :N_reads_sum)
grouped = select(df, [:Label,:tax_name,:N_reads,:damage])

# Compute total reads per sample
total_reads = combine(groupby(df, :Label), :N_reads => sum => :total_N_reads)

# Join total reads to the grouped table
grouped = leftjoin(grouped, total_reads, on = :Label)

# Compute relative abundance
grouped.relative_abundance = grouped.N_reads ./ grouped.total_N_reads

# Pivot to wide format: rows = samples, columns = taxa
taxa_matrix = @pipe unstack(grouped, :Label, :tax_name, :relative_abundance) |> coalesce.(_, 0)
dmg_matrix = @pipe unstack(grouped, :Label, :tax_name, :damage) 

#Step 2: Create sample_data Table
# Read in the DNA concentration CSV
mt = CSV.read("./metadata/HIBO_library_metadata.csv", DataFrame)

# Determine if each Label in df is a sample or control
labels = unique(df.Label)
sample_status = DataFrame(Label = labels)
sample_status.sample_type = ifelse.(in.(sample_status.Label, Ref(x.Label)), "sample", "control")

# Merge with DNA concentrations
sample_data = leftjoin(sample_status, select(mt,["Label","molarity"]), on = :Label)

CSV.write("./deContamination/output/ra_matrix.$rank.csv", taxa_matrix)
CSV.write("./deContamination/output/dmg_matrix.$rank.csv", dmg_matrix)
CSV.write("./deContamination/output/sample_data.$rank.csv", sample_data)






