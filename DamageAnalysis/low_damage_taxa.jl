#!/usr/bin/env julia
using DataFrames
using JLD2
using CSV #CSVFiles
using Pipe
using Formatting

include("./MicrobeProfiling/_get_annotated_df.jl")


pip="amaw"; alg=".lca";add="_ANI92";frank="species";
tag=pip*alg*add

df = load_object("./deContamination/data/$tag.$frank.jld2")

# filter for non-ancient taxa using the pre-calculated thresholds for each sample
stats = load("./DamageAnalysis/output_ANI95/amaw_Bacteria_species_filter_thres_after_ANI95_strict.csv") |> DataFrame 

libs = load("./metadata/HIBO_library_metadata.csv") |> DataFrame |> dropmissing

df_with_thresholds = leftjoin(df, select(stats, [:Label, :lower]), on = :Label)

# Then filter for active microbes and add tax info
filtered_df = @pipe filter(row -> row.damage <= row.lower, df_with_thresholds) |>
              filter(:N_reads=>n->n>1000, _) |> 
              DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank(p, "phylum"))=> :phylum) |>
              DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "class"))=> :class) |>
              DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "family"))=> :family) |>
              DataFrames.transform(_, :tax_path=> ByRow(p->tax_at_rank.(p, "genus"))=> :genus)

# Remove the threshold column 
df_simple = select(filtered_df, [:Label,:Middle_depth, :phylum, :family,:genus, :N_reads])

# write to csv
CSV.write("./DamageAnalysis/output/active_taxa.csv", df_simple)

combine(groupby(filtered_df, :Label), :N_reads=>sum=>:total_reads)
combine(groupby(filtered_df, [:Label]), nrow => :count)
